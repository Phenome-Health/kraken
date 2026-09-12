"""Build the weighted match graph: CURIE pairs with accumulated evidence weight.

Evidence comes from three places (plan §1):

1. **Equivalency lists -> cliques** (capped; star beyond the cap).
2. **Source match predicates** (exact_match/same_as full, close_match low;
   broad_match/narrow_match excluded).
3. **Name similarity** (primary names only; see ``name_norm``).

Each evidence item is tagged with a *source group* so accumulation can
de-correlate the aggregators: within a group weights combine by **max**, across
groups by **sum** (plan §1, and the ACE/RTD correlated-aggregator problem).

The public output is an iterator of ``(a, b, weight)`` with ``a < b`` and
``weight >= tau``. For build scale (~30M CURIEs) the accumulation should move to
the on-disk keyed-sort pattern already used by ``integrate.integrate_edges``;
the in-memory accumulator here is correct and used for smaller runs and tests.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from kraken.entity_resolution.weights import ERWeights

# One piece of evidence for a pair: (a, b, source_group, weight), a < b.
Evidence = tuple[str, str, str, float]
WeightedPair = tuple[str, str, float]


def _ordered(a: str, b: str) -> tuple[str, str]:
    return (a, b) if a < b else (b, a)


def clique_evidence(
    equivalent_ids: Iterable[str],
    source: str,
    weights: ERWeights,
) -> Iterator[Evidence]:
    """Emit evidence edges for one source's equivalency set.

    Up to ``clique_cap`` ids -> full clique (robust: survives to the accumulated
    weight-vs-gamma threshold, all-or-nothing). Beyond the cap -> a star from the
    lexically smallest id (fragile on purpose; big "equivalent" lists are junk).

    The weight depends on the list's SIZE for size-aware sources (see
    ``ERWeights.max_merge_list_size``); a list large enough to carry no weight at all
    emits nothing, rather than a pile of zero-weight edges. Its ids are still seeded as
    nodes elsewhere -- only the equivalence claim is dropped.
    """
    ids = sorted({i for i in equivalent_ids if i})
    if len(ids) < 2:
        return
    group = weights.source_group(source)
    weight = weights.equivalency_weight(source, len(ids))
    if weight <= 0:
        return
    if len(ids) <= weights.clique_cap:
        n = len(ids)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = _ordered(ids[i], ids[j])
                yield (a, b, group, weight)
    else:
        hub = ids[0]
        for other in ids[1:]:
            a, b = _ordered(hub, other)
            yield (a, b, group, weight)


def alias_evidence(
    original: str,
    canonical: str,
    source: str,
    weights: ERWeights,
) -> Evidence | None:
    """Evidence for one ``original -> canonical`` alias asserted by a canonicalized aggregator edge.

    Weighted exactly like a two-id equivalency list from the same source (``ERWeights.alias_weight``)
    -- it is the same kind of claim -- and kept in that source's normal group, so an aggregator
    echoing what the NN clique already says still counts once (max within ``sri_nn_derived``).
    """
    if not original or not canonical or original == canonical:
        return None
    a, b = _ordered(original, canonical)
    return (a, b, weights.source_group(source), weights.alias_weight(source))


def match_predicate_evidence(
    subject: str,
    object_: str,
    predicate: str,
    source: str,
    weights: ERWeights,
    primary_ks: str | None = None,
) -> Evidence | None:
    """Evidence for one source match-predicate edge, or None if it must not
    contribute (hierarchical / non-match predicate / self-loop).

    Grouped per (source, primary knowledge source), so parallel assertions of the same
    pair from different primary KSes sum instead of collapsing to their max (see
    ``ERWeights.predicate_group``)."""
    weight = weights.predicate_weight(predicate)
    if weight is None or not subject or not object_ or subject == object_:
        return None
    a, b = _ordered(subject, object_)
    return (a, b, weights.predicate_group(source, primary_ks), weight)


def name_similarity_evidence(
    pairs: Iterable[tuple[str, str]],
    weights: ERWeights,
) -> Iterator[Evidence]:
    """Wrap name-similarity pairs (from ``name_norm``) as evidence. Name evidence
    is its own source group so it never double-counts with itself."""
    from kraken.entity_resolution.weights import NAME_SIMILARITY_GROUP

    weight = weights.name_similarity_weight
    for a, b in pairs:
        if a == b:
            continue
        a, b = _ordered(a, b)
        yield (a, b, NAME_SIMILARITY_GROUP, weight)


def accumulate(evidence: Iterable[Evidence], weights: ERWeights) -> dict[tuple[str, str], float]:
    """Accumulate evidence into per-pair total weight.

    Within a source group: **max**. Across groups: **sum**. This treats
    correlated sources (the SRI-NN-derived aggregators) as a single source while
    letting independent sources reinforce each other.
    """
    # pair -> {group -> max weight seen}
    per_pair: dict[tuple[str, str], dict[str, float]] = {}
    for a, b, group, weight in evidence:
        key = (a, b)
        groups = per_pair.get(key)
        if groups is None:
            per_pair[key] = {group: weight}
        else:
            prev = groups.get(group)
            if prev is None or weight > prev:
                groups[group] = weight
    return {pair: sum(groups.values()) for pair, groups in per_pair.items()}


def filter_by_tau(
    totals: dict[tuple[str, str], float],
    weights: ERWeights,
) -> Iterator[WeightedPair]:
    """Yield ``(a, b, weight)`` for pairs meeting the tau pre-filter, sorted for
    determinism."""
    for (a, b), weight in sorted(totals.items()):
        if weight >= weights.tau:
            yield (a, b, weight)


@dataclass
class MatchGraph:
    """Convenience wrapper: accumulate evidence and expose the tau-filtered graph."""

    weights: ERWeights

    def build(self, evidence: Iterable[Evidence]) -> list[WeightedPair]:
        totals = accumulate(evidence, self.weights)
        return list(filter_by_tau(totals, self.weights))
