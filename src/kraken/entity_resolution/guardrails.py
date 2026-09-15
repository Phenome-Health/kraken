"""Guardrails: a safety net checked on clustering output (plan §3).

Every guardrail is **hereditary** (closed under subsets), so "split until valid"
is well-founded: singletons are always valid, so the recursion terminates.

Guardrails implemented:

* **One Biolink branch** — cluster valid iff the intersection of ``branches``
  across nodes is non-empty (see ``families.BranchFamilies``).
* **One taxon** — at most one distinct taxon among nodes that have one (untaxoned
  nodes are wildcards). Taxa are species-normalized at harmonization.
* **One id per structural prefix** — RefMet / LIPID MAPS ids each denote one
  structural entity by construction. Candidate prefixes (HGNC, NCBIGene) are
  *instrumented but not enforced* by default; promote once the histogram is
  clean.

Repair strategy for a violating cluster: an optional injected ``splitter`` may try
to split it better; otherwise (or if it can't) fall back to a deterministic greedy
valid partition that regrows the cluster along its strongest edges, keeping only
merges that stay valid (repairs but cannot discover). A one-id violation is always
repaired, however many ids are involved: k ids of a one-entity-per-id prefix in one
cluster means k entities were merged, so k clusters is the right answer, and a large
k is the strongest sign of it -- not a reason to stop. Large repairs are logged,
because they point at an upstream conflation worth finding. (This used to be capped
at 3, leaving bigger violations intact; that shipped the worst conflations
unrepaired, e.g. a single 2.1.1 cluster holding 247 RefMet and 247 LIPID MAPS ids.)

Two documented blind spots (log, don't solve): a node violating a rule *by
itself* (evaluated cross-node only), and protein vs. cleavage products (Biolink
types them inconsistently).
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass

from kraken.entity_resolution.families import ALL_FAMILIES, BranchFamilies
from kraken.utils.constants import SMILES_PREFIX

# Splitter injected by resolve/clustering: given member ids, return sub-clusters.
Splitter = Callable[[list[str]], list[list[str]]]

# One-id-per-cluster prefixes ENFORCED by default:
#   * RM / LM — RefMet and LIPID MAPS, each one structural entity by construction
#     (verified prefixes; NOT "REFMET"/"LIPIDMAPS").
#   * MONDO — curated one-id-per-disease authority. Watch the histogram:
#     obsoleted-and-replaced MONDO terms can legitimately co-occur, so this can
#     cause false splits. (2.1.1 had exactly one cluster with >1 MONDO id.) If the
#     histogram shows false splits, demote MONDO to a candidate prefix rather than
#     weakening the repair for RM/LM, which have no legitimate co-occurrence.
#   * CAID — ClinGen Allele Registry: one canonical allele per id, so two in a cluster
#     means two alleles were merged (e.g. by a shared rsid, which names the position
#     they sit at rather than either allele).
#   * SMILES — one structure per id once biomapper2 has canonicalized it. Two in a cluster are two structures:
#     in practice a compound and its salt, or two stereoisomers (of 846 lipids lipidmaps and translator share an
#     InChIKey with, the 60 whose SMILES differ are exactly that). Drug/chemical conflation is off, so they stay
#     apart.
DEFAULT_ENFORCED_PREFIXES: frozenset[str] = frozenset({"RM", "LM", "MONDO", "CAID", SMILES_PREFIX})

# Candidate one-id-per-cluster prefixes: watched (instrumented) but NOT enforced.
# HGNC only — it is a curated one-id-per-human-gene nomenclature, so >1 usually
# means two genes were merged. Promote to enforced only once the ids-per-cluster
# histogram is clean (paralogs that share a protein can legitimately co-occur
# under gene/protein conflation; HGNC also has withdrawn/replaced ids). NCBIGene
# is deliberately NOT a candidate: it spans all species, carries
# obsoleted-and-replaced ids, and correctly co-occurs when paralogs share a
# protein under gene/protein conflation.
DEFAULT_CANDIDATE_PREFIXES: frozenset[str] = frozenset({"HGNC"})


@dataclass(frozen=True)
class NodeInfo:
    """Per-CURIE facts the guardrails need. ``branches`` is the node's precomputed
    branch-family set (resolved from its category upstream, in build), so every
    guardrail check operates on families directly — no repeated category->family
    conversion during the split-until-valid recursion. Default is the wildcard
    ``ALL_FAMILIES`` (an untyped node never constrains a merge)."""

    curie: str
    branches: frozenset[str] = ALL_FAMILIES
    taxon: str | None = None

    @property
    def prefix(self) -> str:
        return self.curie.split(":", 1)[0]


@dataclass
class GuardrailConfig:
    enforced_prefixes: frozenset[str] = DEFAULT_ENFORCED_PREFIXES
    candidate_prefixes: frozenset[str] = DEFAULT_CANDIDATE_PREFIXES
    # A one-id repair involving more ids of one prefix than this is logged (it is
    # still carried out) -- a large one points at an upstream conflation.
    one_id_repair_log_threshold: int = 3
    # Log clusters at or above this size (no hard maximum — no defensible one).
    oversized_cluster_log_threshold: int = 100


NodeInfoMap = Mapping[str, NodeInfo]


def branch_valid(members: Iterable[str], info: NodeInfoMap) -> bool:
    node_branches = [info[m].branches if m in info else ALL_FAMILIES for m in members]
    return BranchFamilies.cluster_branches(node_branches) is not None


def taxon_valid(members: Iterable[str], info: NodeInfoMap) -> bool:
    taxa = {info[m].taxon for m in members if m in info and info[m].taxon}
    return len(taxa) <= 1


def one_id_valid(members: Iterable[str], enforced_prefixes: frozenset[str]) -> bool:
    by_prefix: Counter[str] = Counter()
    for m in members:
        p = m.split(":", 1)[0]
        if p in enforced_prefixes:
            by_prefix[p] += 1
    return all(count <= 1 for count in by_prefix.values())


def cluster_violations(
    members: list[str],
    info: NodeInfoMap,
    config: GuardrailConfig,
) -> list[str]:
    """Return the list of guardrail names a cluster violates (empty = valid)."""
    if len(members) <= 1:
        return []
    violations = []
    if not branch_valid(members, info):
        violations.append("branch")
    if not taxon_valid(members, info):
        violations.append("taxon")
    if not one_id_valid(members, config.enforced_prefixes):
        violations.append("one_id")
    return violations


def greedy_valid_partition(
    members: list[str],
    info: NodeInfoMap,
    config: GuardrailConfig,
    adjacency: Mapping[str, Mapping[str, float]],
) -> list[list[str]]:
    """Deterministically partition members into guardrail-valid groups by growing them along the match graph,
    strongest edges first.

    Every member starts as its own group. Edges between members are taken in order of decreasing weight (ties by
    id), and each joins its two endpoints' groups unless the merged group would break a guardrail -- two taxa, two
    ids of one enforced prefix, or no family in common. So the strongest evidence in the cluster is kept, a group
    only ever grows along an edge (every group is connected), and a member none of whose edges can be kept stays
    on its own rather than being put somewhere it has no evidence for. Every group is valid by construction.

    Growing in id order instead -- the previous approach -- could place a member before any of its neighbours, drop
    it into the earliest group that allowed it, and so use up that group's one MONDO slot ahead of the MONDO id
    the group was actually connected to (MONDO:1 -0.5- UMLS:C1 -0.5- MONDO:2 -1.0- DOID:7 came out as
    [DOID:7, MONDO:1, UMLS:C1] + [MONDO:2], cutting the one strong link).

    Near-linear: one sort of the cluster's edges, then union-find over groups that carry their guardrail state
    (shared families, taxon, enforced prefixes present), so each merge check is O(1) and state merges small into
    large.
    """
    enforced = config.enforced_prefixes
    member_set = set(members)
    parent = {m: m for m in members}
    size = dict.fromkeys(members, 1)
    group_branches: dict[str, frozenset[str] | None] = {}  # None = unconstrained (only wildcards so far)
    group_taxon: dict[str, str | None] = {}
    group_prefixes: dict[str, set[str]] = {}
    for m in members:
        ni = info.get(m)
        branches = ni.branches if ni is not None else ALL_FAMILIES
        group_branches[m] = None if (branches is ALL_FAMILIES or branches == ALL_FAMILIES) else branches
        group_taxon[m] = ni.taxon if ni is not None else None
        prefix = m.split(":", 1)[0]
        group_prefixes[m] = {prefix} if prefix in enforced else set()

    def find(m: str) -> str:
        while parent[m] != m:
            parent[m] = parent[parent[m]]
            m = parent[m]
        return m

    # Each member pair once, whichever side's adjacency lists it.
    pair_weights: dict[tuple[str, str], float] = {}
    for a in members:
        for b, weight in adjacency.get(a, {}).items():
            if b in member_set and b != a:
                pair = (a, b) if a < b else (b, a)
                pair_weights[pair] = max(weight, pair_weights.get(pair, weight))
    edges = sorted((-weight, a, b) for (a, b), weight in pair_weights.items())
    for _neg_weight, a, b in edges:
        root_a, root_b = find(a), find(b)
        if root_a == root_b:
            continue
        taxon_a, taxon_b = group_taxon[root_a], group_taxon[root_b]
        if taxon_a is not None and taxon_b is not None and taxon_a != taxon_b:
            continue
        if not group_prefixes[root_a].isdisjoint(group_prefixes[root_b]):
            continue
        branches_a, branches_b = group_branches[root_a], group_branches[root_b]
        if branches_a is None:
            merged_branches = branches_b
        elif branches_b is None:
            merged_branches = branches_a
        else:
            merged_branches = branches_a & branches_b
        if merged_branches is not None and not merged_branches:
            continue
        if size[root_a] < size[root_b]:
            root_a, root_b = root_b, root_a
        parent[root_b] = root_a
        size[root_a] += size[root_b]
        group_branches[root_a] = merged_branches
        group_taxon[root_a] = taxon_a if taxon_a is not None else taxon_b
        group_prefixes[root_a] |= group_prefixes[root_b]

    groups: dict[str, list[str]] = defaultdict(list)
    for m in members:
        groups[find(m)].append(m)
    return sorted((sorted(group) for group in groups.values()), key=lambda group: group[0])


def enforce_cluster(
    members: list[str],
    info: NodeInfoMap,
    config: GuardrailConfig,
    *,
    adjacency: Mapping[str, Mapping[str, float]],
    splitter: Splitter | None = None,
) -> list[list[str]]:
    """Split a cluster until every part is guardrail-valid.

    Tries ``splitter`` first if provided; otherwise (or if it fails to reduce the
    cluster) falls back to ``greedy_valid_partition``. Terminates because both
    fallbacks strictly shrink clusters and singletons are always valid.
    """
    members = sorted(members)
    violations = cluster_violations(members, info, config)
    if not violations:
        return [members]

    # Surface large one-id repairs: they mean many distinct entities were merged upstream.
    if "one_id" in violations:
        n_ids = _max_offending_id_count(members, config.enforced_prefixes)
        if n_ids > config.one_id_repair_log_threshold:
            logging.warning(
                "one_id violation with %d ids of a single prefix; splitting the cluster (%d members) -- "
                "a repair this large usually means an upstream conflation (members=%s...)",
                n_ids,
                len(members),
                members[:6],
            )

    sub: list[list[str]] | None = None
    if splitter is not None:
        candidate = splitter(members)
        if len(candidate) > 1:
            sub = candidate
    if sub is None:
        sub = greedy_valid_partition(members, info, config, adjacency)
    if len(sub) <= 1:
        # nothing split it (e.g. an unsplittable single-branch blob) -> stop
        logging.warning("could not split violating cluster %s (violations=%s)", members[:6], violations)
        return [members]

    result: list[list[str]] = []
    for part in sub:
        if part == members:  # no progress; avoid infinite recursion
            result.append(part)
        else:
            result.extend(enforce_cluster(part, info, config, splitter=splitter, adjacency=adjacency))
    return result


def _max_offending_id_count(members: Iterable[str], enforced_prefixes: frozenset[str]) -> int:
    by_prefix: Counter[str] = Counter()
    for m in members:
        p = m.split(":", 1)[0]
        if p in enforced_prefixes:
            by_prefix[p] += 1
    return max(by_prefix.values(), default=0)


# ---- instrumentation (plan: emit before promoting a candidate rule) ----


def ids_per_cluster_histogram(clusters: Iterable[Iterable[str]]) -> dict[str, dict[int, int]]:
    """Per-prefix histogram of {ids-of-that-prefix-in-a-cluster: number of such
    clusters}. Feeds the decision to promote a candidate one-id rule."""
    hist: dict[str, Counter[int]] = defaultdict(Counter)
    for cluster in clusters:
        by_prefix: Counter[str] = Counter()
        for m in cluster:
            by_prefix[m.split(":", 1)[0]] += 1
        for prefix, count in by_prefix.items():
            hist[prefix][count] += 1
    return {prefix: dict(counter) for prefix, counter in hist.items()}


def log_oversized_clusters(clusters: Iterable[Iterable[str]], config: GuardrailConfig) -> list[list[str]]:
    """Log (don't reject) clusters at/above the size threshold. Returns them."""
    materialized = [sorted(c) for c in clusters]
    oversized = [c for c in materialized if len(c) >= config.oversized_cluster_log_threshold]
    for cluster in oversized:
        logging.warning("oversized cluster (%d members): %s...", len(cluster), cluster[:8])
    return oversized
