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

Repair strategy for a violating cluster: first try to split it better (raise
gamma and re-Leiden on the induced subgraph — injected as ``splitter``); if that
doesn't split it, fall back to a deterministic greedy valid partition that
respects edge connectivity (seeded label propagation's role — repairs but cannot
discover). The one-id repair is capped: forcing k clusters for k ids is fine at 2
but absurd at 6, so beyond the cap we log and leave the cluster intact.

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
#     cause false splits; the capped repair keeps it bounded.
DEFAULT_ENFORCED_PREFIXES: frozenset[str] = frozenset({"RM", "LM", "MONDO"})

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
    """Per-CURIE facts the guardrails need."""

    curie: str
    categories: tuple[str, ...] = ()
    taxon: str | None = None

    @property
    def prefix(self) -> str:
        return self.curie.split(":", 1)[0]


@dataclass
class GuardrailConfig:
    enforced_prefixes: frozenset[str] = DEFAULT_ENFORCED_PREFIXES
    candidate_prefixes: frozenset[str] = DEFAULT_CANDIDATE_PREFIXES
    # Cap the one-id repair: don't force more than this many sub-clusters to
    # satisfy a one-id rule (log oversized instead).
    one_id_repair_cap: int = 3
    # Log clusters at or above this size (no hard maximum — no defensible one).
    oversized_cluster_log_threshold: int = 100


NodeInfoMap = Mapping[str, NodeInfo]


def _branches_of(members: Iterable[str], info: NodeInfoMap, families: BranchFamilies) -> list[frozenset[str]]:
    return [families.branches(info[m].categories if m in info else None) for m in members]


def branch_valid(members: Iterable[str], info: NodeInfoMap, families: BranchFamilies) -> bool:
    return BranchFamilies.cluster_branches(_branches_of(members, info, families)) is not None


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
    families: BranchFamilies,
    config: GuardrailConfig,
) -> list[str]:
    """Return the list of guardrail names a cluster violates (empty = valid)."""
    if len(members) <= 1:
        return []
    violations = []
    if not branch_valid(members, info, families):
        violations.append("branch")
    if not taxon_valid(members, info):
        violations.append("taxon")
    if not one_id_valid(members, config.enforced_prefixes):
        violations.append("one_id")
    return violations


def _node_valid_in_group(
    node: str,
    group: list[str],
    info: NodeInfoMap,
    families: BranchFamilies,
    config: GuardrailConfig,
) -> bool:
    """Would adding ``node`` to ``group`` keep every guardrail satisfied?"""
    candidate = group + [node]
    return not cluster_violations(candidate, info, families, config)


def greedy_valid_partition(
    members: list[str],
    info: NodeInfoMap,
    families: BranchFamilies,
    config: GuardrailConfig,
    adjacency: Mapping[str, Mapping[str, float]] | None = None,
) -> list[list[str]]:
    """Deterministically partition members into guardrail-valid groups, placing
    each node into the connected-most valid group (new group if none fits).

    Constraining (non-wildcard) nodes are processed first so they seed distinct
    groups; wildcard nodes then attach by connectivity (plan: "placed by
    connectivity when a blob splits").
    """
    adjacency = adjacency or {}

    def is_wildcard(m: str) -> bool:
        if m not in info:
            return True
        ni = info[m]
        branches = families.branches(ni.categories)
        return (
            (branches is ALL_FAMILIES or branches == ALL_FAMILIES)
            and ni.taxon is None
            and ni.prefix not in (config.enforced_prefixes | config.candidate_prefixes)
        )

    ordered = sorted(members, key=lambda m: (is_wildcard(m), m))
    groups: list[list[str]] = []
    for node in ordered:
        best_idx = -1
        best_score = None
        for idx, group in enumerate(groups):
            if not _node_valid_in_group(node, group, info, families, config):
                continue
            score = sum(adjacency.get(node, {}).get(other, 0.0) for other in group)
            # prefer higher connectivity; tie-break to earliest group (deterministic)
            if best_score is None or score > best_score:
                best_score = score
                best_idx = idx
        if best_idx >= 0:
            groups[best_idx].append(node)
        else:
            groups.append([node])
    return [sorted(g) for g in groups]


def enforce_cluster(
    members: list[str],
    info: NodeInfoMap,
    families: BranchFamilies,
    config: GuardrailConfig,
    *,
    splitter: Splitter | None = None,
    adjacency: Mapping[str, Mapping[str, float]] | None = None,
) -> list[list[str]]:
    """Split a cluster until every part is guardrail-valid.

    Tries ``splitter`` (raise gamma + re-Leiden) first; if it fails to reduce the
    cluster, falls back to ``greedy_valid_partition``. Terminates because both
    fallbacks strictly shrink clusters and singletons are always valid.
    """
    members = sorted(members)
    violations = cluster_violations(members, info, families, config)
    if not violations:
        return [members]

    # Guard against the one-id repair blowing up (plan: cap the repair).
    if violations == ["one_id"]:
        n_ids = _max_offending_id_count(members, config.enforced_prefixes)
        if n_ids > config.one_id_repair_cap:
            logging.warning(
                "one_id violation with %d ids of a single prefix exceeds repair cap %d; "
                "leaving cluster intact and logging (members=%s...)",
                n_ids,
                config.one_id_repair_cap,
                members[:6],
            )
            return [members]

    sub: list[list[str]] | None = None
    if splitter is not None:
        candidate = splitter(members)
        if len(candidate) > 1:
            sub = candidate
    if sub is None:
        sub = greedy_valid_partition(members, info, families, config, adjacency)
    if len(sub) <= 1:
        # nothing split it (e.g. an unsplittable single-branch blob) -> stop
        logging.warning("could not split violating cluster %s (violations=%s)", members[:6], violations)
        return [members]

    result: list[list[str]] = []
    for part in sub:
        if part == members:  # no progress; avoid infinite recursion
            result.append(part)
        else:
            result.extend(enforce_cluster(part, info, families, config, splitter=splitter, adjacency=adjacency))
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
