"""Tests for match-graph construction, accumulation, and source-group de-correlation."""

import pytest

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
    # sri_nn_derived source group -> combine by MAX, not sum.
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


# --------------------------------------------------------------------------------------
# Size-aware equivalency weights
# --------------------------------------------------------------------------------------


def test_aggregator_weight_steps_down_with_list_size():
    """kg2: merge-strength up to 24 ids, corroboration-only to 60, nothing beyond."""
    w = ERWeights()
    assert w.equivalency_weight("kg2", 2) >= w.tau
    assert w.equivalency_weight("kg2", 24) >= w.tau
    assert 0 < w.equivalency_weight("kg2", 25) < w.tau
    assert 0 < w.equivalency_weight("kg2", 60) < w.tau
    assert w.equivalency_weight("kg2", 61) == 0


def test_robokop_merges_larger_lists_than_kg2():
    """robokop's lists degrade more gracefully (87% vs kg2's 73% at 21-30), so its threshold is higher."""
    w = ERWeights()
    assert w.equivalency_weight("robokop", 30) >= w.tau
    assert w.equivalency_weight("kg2", 30) < w.tau
    assert w.equivalency_weight("robokop", 31) < w.tau


@pytest.mark.parametrize("source", ["nn", "ncbigene", "umls", "refmet"])
def test_sources_not_listed_as_size_aware_keep_a_flat_weight(source):
    """The size curve was measured on aggregators only. The normalizer's cliques in particular are
    clean and legitimately large (gene/protein ~30-40); capping them would break up the backbone."""
    w = ERWeights()
    assert w.equivalency_weight(source, 2) == w.equivalency_weight(source, 40) == w.equivalency_weight(source, 500)


def test_alias_is_weighted_as_a_two_id_list():
    """An alias is the same kind of claim as a list entry, so the two can never drift apart."""
    w = ERWeights()
    for source in ["kg2", "robokop", "translator-kg-open"]:
        assert w.alias_weight(source) == w.equivalency_weight(source, 2)
        assert w.alias_weight(source) >= w.tau


def test_corroborate_weight_must_stay_below_tau():
    """Otherwise mid-size lists would merge on their own, defeating the size-awareness."""
    with pytest.raises(ValueError, match="corroborate_weight"):
        ERWeights(corroborate_weight=0.3, tau=0.3)


def test_list_too_large_to_carry_weight_emits_no_evidence():
    """Not a pile of zero-weight edges -- nothing at all (its ids are still seeded as nodes elsewhere)."""
    w = ERWeights()
    ids = [f"C:{i}" for i in range(61)]
    assert list(clique_evidence(ids, "kg2", w)) == []
    assert list(clique_evidence(ids[:60], "kg2", w))  # one fewer still corroborates


def test_every_clique_edge_carries_the_whole_list_size_weight():
    """A 25-id kg2 list is judged as a 25-id list on every one of its edges, not per pair."""
    w = ERWeights()
    ev = list(clique_evidence([f"C:{i}" for i in range(25)], "kg2", w))
    assert {wt for _a, _b, _g, wt in ev} == {w.equivalency_weight("kg2", 25)}


# --------------------------------------------------------------------------------------
# Parallel match predicates from different primary knowledge sources
# --------------------------------------------------------------------------------------


def test_parallel_close_matches_from_different_primary_kses_sum_to_a_merge():
    """Three independent KSes each asserting close_match on one pair are three claims, not one."""
    w = ERWeights()
    ev = [
        match_predicate_evidence("A:1", "B:1", "biolink:close_match", "kg2", w, primary_ks=ks)
        for ks in ("infores:mesh", "infores:go", "infores:chv-umls")
    ]
    total = accumulate(ev, w)[("A:1", "B:1")]
    assert total == pytest.approx(3 * w.close_match_weight)
    assert total >= w.tau


def test_parallel_close_matches_from_the_same_primary_ks_do_not_sum():
    """Repeating the same KS's claim is not corroboration."""
    w = ERWeights()
    ev = [
        match_predicate_evidence("A:1", "B:1", "biolink:close_match", "kg2", w, primary_ks="infores:mesh")
        for _ in range(3)
    ]
    assert accumulate(ev, w)[("A:1", "B:1")] == w.close_match_weight


def test_close_match_without_a_primary_ks_falls_back_to_the_source_group():
    """Unattributed claims from one source must not count as independent of each other."""
    w = ERWeights()
    ev = [match_predicate_evidence("A:1", "B:1", "biolink:close_match", "kg2", w) for _ in range(3)]
    assert accumulate(ev, w)[("A:1", "B:1")] == w.close_match_weight
    assert ev[0][2] == w.source_group("kg2")
