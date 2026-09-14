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

Repair strategy for a violating cluster: an optional injected ``splitter`` may
try to split it better; otherwise (or if it can't) fall back to a deterministic
greedy valid partition that respects edge connectivity (repairs but cannot
discover). A one-id violation is always repaired, however many ids are involved:
k ids of a one-entity-per-id prefix in one cluster means k entities were merged,
so k clusters is the right answer, and a large k is the strongest sign of it --
not a reason to stop. Large repairs are logged, because they point at an upstream
conflation worth finding. (This used to be capped at 3, leaving bigger violations
intact; that shipped the worst conflations unrepaired, e.g. a single 2.1.1 cluster
holding 247 RefMet and 247 LIPID MAPS ids.)

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
DEFAULT_ENFORCED_PREFIXES: frozenset[str] = frozenset({"RM", "LM", "MONDO", "CAID"})

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
    adjacency: Mapping[str, Mapping[str, float]] | None = None,
) -> list[list[str]]:
    """Deterministically partition members into guardrail-valid groups, placing
    each node into the connected-most valid group (new group if none fits; ties go
    to the earliest group).

    Constraining (non-wildcard) nodes are processed first so they seed distinct
    groups; wildcard nodes then attach by connectivity (plan: "placed by
    connectivity when a blob splits").

    Near-linear in the cluster's size. Each group keeps its guardrail state (branch
    intersection, taxon, enforced prefixes present) so a node is checked against a
    group in O(1), connectivity is scored only over the node's own neighbours, and
    the earliest compatible group is found with a per-node-kind pointer that only
    moves forward: a group only becomes more constrained as members join, so once it
    rejects a kind of node it rejects that kind for good. (It used to rebuild and
    re-validate every candidate group for every node -- O(n^2), which became
    reachable once large one-id repairs stopped being capped.)
    """
    adjacency = adjacency or {}
    enforced = config.enforced_prefixes

    def is_wildcard(m: str) -> bool:
        if m not in info:
            return True
        ni = info[m]
        return (
            (ni.branches is ALL_FAMILIES or ni.branches == ALL_FAMILIES)
            and ni.taxon is None
            and ni.prefix not in (enforced | config.candidate_prefixes)
        )

    def kind_of(m: str) -> tuple[frozenset[str] | None, str | None, str | None]:
        """(branches, or None for a wildcard; taxon; enforced prefix, or None)."""
        ni = info.get(m)
        branches = ni.branches if ni is not None else ALL_FAMILIES
        prefix = m.split(":", 1)[0]
        return (
            None if (branches is ALL_FAMILIES or branches == ALL_FAMILIES) else branches,
            ni.taxon if ni is not None else None,
            prefix if prefix in enforced else None,
        )

    groups: list[list[str]] = []
    group_branches: list[frozenset[str] | None] = []  # running intersection; None = no constraint yet
    group_taxon: list[str | None] = []
    group_prefixes: list[set[str]] = []
    group_of: dict[str, int] = {}
    earliest_open: dict[tuple, int] = {}  # node kind -> first group index that may still accept it

    def accepts(idx: int, branches: frozenset[str] | None, taxon: str | None, prefix: str | None) -> bool:
        if prefix is not None and prefix in group_prefixes[idx]:
            return False
        if taxon is not None and group_taxon[idx] is not None and group_taxon[idx] != taxon:
            return False
        current = group_branches[idx]
        if current is not None and not current:  # a member with no branches: nothing may join it
            return False
        if branches is not None:
            if not branches:  # a node with no branches can't join anything
                return False
            if current is not None and not (current & branches):
                return False
        return True

    for node in sorted(members, key=lambda m: (is_wildcard(m), m)):
        kind = kind_of(node)
        branches, taxon, prefix = kind

        # connectivity scores, over the node's neighbours only
        scores: dict[int, float] = defaultdict(float)
        for neighbour, weight in adjacency.get(node, {}).items():
            idx = group_of.get(neighbour)
            if idx is not None:
                scores[idx] += weight

        # the earliest group that accepts this kind of node (pointer only moves forward)
        idx = earliest_open.get(kind, 0)
        while idx < len(groups) and not accepts(idx, *kind):
            idx += 1
        earliest_open[kind] = idx

        best_idx, best_score = -1, None
        if idx < len(groups):
            best_idx, best_score = idx, scores.get(idx, 0.0)
        for candidate, score in scores.items():
            if score > (best_score if best_score is not None else float("-inf")) or (
                score == best_score and candidate < best_idx
            ):
                if accepts(candidate, *kind):
                    best_idx, best_score = candidate, score

        if best_idx < 0:
            best_idx = len(groups)
            groups.append([])
            group_branches.append(None)
            group_taxon.append(None)
            group_prefixes.append(set())
        groups[best_idx].append(node)
        group_of[node] = best_idx
        if branches is not None:
            current = group_branches[best_idx]
            group_branches[best_idx] = branches if current is None else current & branches
        if taxon is not None:
            group_taxon[best_idx] = taxon
        if prefix is not None:
            group_prefixes[best_idx].add(prefix)

    return [sorted(g) for g in groups]


def enforce_cluster(
    members: list[str],
    info: NodeInfoMap,
    config: GuardrailConfig,
    *,
    splitter: Splitter | None = None,
    adjacency: Mapping[str, Mapping[str, float]] | None = None,
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
