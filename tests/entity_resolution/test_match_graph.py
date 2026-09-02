"""Tests for match-graph construction, accumulation, and de-correlation."""

from kraken.entity_resolution.match_graph import (
    accumulate,
    clique_evidence,
    match_predicate_evidence,
    name_similarity_evidence,
)
from kraken.entity_resolution.weights import ERWeights


def test_clique_evidence_full_below_cap():
    w = ERWeights(clique_cap=20)
    ev = list(clique_evidence(["C:3", "C:1", "C:2"], "ncbigene", w))
    pairs = {(a, b) for a, b, _g, _wt in ev}
    assert pairs == {("C:1", "C:2"), ("C:1", "C:3"), ("C:2", "C:3")}
    # sorted endpoints
    assert all(a < b for a, b, _g, _wt in ev)


def test_clique_evidence_star_above_cap():
    w = ERWeights(clique_cap=3)
    ids = [f"C:{i}" for i in range(10)]
    ev = list(clique_evidence(ids, "ncbigene", w))
    # star from the lexically smallest id -> n-1 edges, all touching the hub
    assert len(ev) == 9
    hub = min(ids)
    assert all(hub in (a, b) for a, b, _g, _wt in ev)


def test_excluded_and_close_predicates():
    w = ERWeights()
    assert match_predicate_evidence("A:1", "B:1", "biolink:broad_match", "kg2", w) is None
    assert match_predicate_evidence("A:1", "B:1", "biolink:narrow_match", "kg2", w) is None
    assert match_predicate_evidence("A:1", "A:1", "biolink:exact_match", "kg2", w) is None  # self loop
    close = match_predicate_evidence("B:1", "A:1", "biolink:close_match", "kg2", w)
    assert close is not None
    a, b, _g, wt = close
    assert (a, b) == ("A:1", "B:1")  # reordered
    assert wt == w.close_match_weight


def test_accumulate_sums_independent_sources():
    w = ERWeights()
    # Two independent curated sources assert the same pair -> weights add.
    ev = list(clique_evidence(["A:1", "B:1"], "ncbigene", w)) + list(clique_evidence(["A:1", "B:1"], "refmet", w))
    totals = accumulate(ev, w)
    assert totals[("A:1", "B:1")] == w.equivalency_weight("ncbigene") + w.equivalency_weight("refmet")


def test_accumulate_decorrelates_aggregators():
    w = ERWeights()
    # KG2, ROBOKOP, Translator all assert the same pair. They share the
    # sri_nn_derived correlation group -> combine by MAX, not sum.
    ev = []
    for src in ["kg2", "robokop", "translator-kg-open"]:
        ev += list(clique_evidence(["A:1", "B:1"], src, w))
    totals = accumulate(ev, w)
    expected = max(w.equivalency_weight(s) for s in ["kg2", "robokop", "translator-kg-open"])
    assert totals[("A:1", "B:1")] == expected


def test_accumulate_mixes_max_within_sum_across():
    w = ERWeights()
    ev = list(clique_evidence(["A:1", "B:1"], "kg2", w))  # aggregator group
    ev += list(clique_evidence(["A:1", "B:1"], "robokop", w))  # same group -> max
    ev += list(clique_evidence(["A:1", "B:1"], "ncbigene", w))  # independent -> sum
    totals = accumulate(ev, w)
    agg = max(w.equivalency_weight("kg2"), w.equivalency_weight("robokop"))
    assert totals[("A:1", "B:1")] == agg + w.equivalency_weight("ncbigene")


def test_name_similarity_own_group():
    w = ERWeights()
    ev = list(name_similarity_evidence([("A:1", "B:1")], w))
    assert len(ev) == 1
    _a, _b, group, wt = ev[0]
    assert wt == w.name_similarity_weight
    # name evidence in its own group doesn't sum against equivalency of same pair
    combined = accumulate(list(clique_evidence(["A:1", "B:1"], "kg2", w)) + ev, w)
    assert combined[("A:1", "B:1")] == w.equivalency_weight("kg2") + w.name_similarity_weight
