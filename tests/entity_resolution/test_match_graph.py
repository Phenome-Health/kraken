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
    # babel_derived source group -> combine by MAX, not sum.
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
# Prefix-capped equivalency lists
# --------------------------------------------------------------------------------------


def _by_pair(evidence):
    return {(a, b): wt for a, b, _g, wt in evidence}


def _pair(a: str, b: str) -> tuple[str, str]:
    """Evidence pairs are emitted in sorted order."""
    return (a, b) if a < b else (b, a)


def test_bulk_prefix_ids_get_only_a_weak_link_to_the_listing_node():
    """TP53's shape: a few good ids plus hundreds from one prefix. The few merge; the bulk only links to TP53."""
    w = ERWeights()
    reactome = [f"REACT:R-HSA-{i}" for i in range(w.max_ids_per_prefix["kg2"] + 1)]
    ids = ["NCBIGene:7157", "HGNC:11998", "NCIT:C17359", *reactome]
    ev = _by_pair(clique_evidence(ids, "kg2", w, head="NCBIGene:7157"))
    for core in ("HGNC:11998", "NCIT:C17359"):  # core: full weight, to the listing node (a star, not a clique)
        assert ev[_pair("NCBIGene:7157", core)] == w.equivalency_weight("kg2")
    assert _pair("HGNC:11998", "NCIT:C17359") not in ev
    for r in reactome:
        assert ev[("NCBIGene:7157", r)] == w.bulk_prefix_weight < w.tau  # bulk: one weak edge to the head
        assert ("HGNC:11998", r) not in ev and ("NCIT:C17359", r) not in ev
    assert not any(a.startswith("REACT:") and b.startswith("REACT:") for a, b in ev)  # never bulk-to-bulk


def test_an_aggregator_list_is_a_star_from_the_listing_node_however_small():
    """An aggregator's list is ONE claim about ONE node. As a clique it is a dense community that label
    propagation prefers over the curated cliques inside it -- which is how kg2's heart list split the heart in
    two, and how its metformin list fused metformin with its hydrochloride."""
    w = ERWeights()
    ids = [f"P{i}:{i}" for i in range(26)]
    ev = _by_pair(clique_evidence(ids, "kg2", w, head="P0:0"))
    assert len(ev) == 25  # a star, not 26*25/2
    assert all("P0:0" in pair for pair in ev)
    assert set(ev.values()) == {w.equivalency_weight("kg2")}


def test_a_curated_source_list_stays_a_full_clique():
    """Babel's cliques and the native curated lists are every-member-vouches-for-every-other, so they stay
    cliques (up to clique_cap) -- that mutual support is what an aggregator list must not imitate."""
    w = ERWeights()
    ids = [f"P{i}:{i}" for i in range(26)]
    ev = list(clique_evidence(ids, "babel", w, head="P0:0"))
    assert len(ev) == 26 * 25 // 2
    assert {wt for _a, _b, _g, wt in ev} == {w.equivalency_weight("babel")}


def test_up_to_the_cap_is_not_bulk():
    w = ERWeights()
    cap = w.max_ids_per_prefix["kg2"]
    ids = ["HEAD:1", *[f"X:{i}" for i in range(cap)]]
    assert {wt for _a, _b, _g, wt in clique_evidence(ids, "kg2", w, head="HEAD:1")} == {w.equivalency_weight("kg2")}


@pytest.mark.parametrize("source", ["babel", "ncbigene", "umls", "refmet"])
def test_uncapped_sources_keep_every_id_at_full_weight(source):
    """Babel's cliques legitimately hold many ids of one prefix (a gene's protein isoforms)."""
    w = ERWeights(clique_cap=1000)
    ids = ["HEAD:1", *[f"ENSEMBL:ENSP{i}" for i in range(50)]]
    assert {wt for _a, _b, _g, wt in clique_evidence(ids, source, w)} == {w.equivalency_weight(source)}


def test_alias_takes_the_sources_full_equivalency_weight():
    w = ERWeights()
    for source in ["kg2", "robokop", "translator-kg-open"]:
        assert w.alias_weight(source) == w.equivalency_weight(source) >= w.tau


def test_bulk_prefix_weight_must_stay_below_tau():
    with pytest.raises(ValueError, match="bulk_prefix_weight"):
        ERWeights(bulk_prefix_weight=0.3, tau=0.3)


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


def test_star_hub_is_the_head_not_the_lexically_smallest_id():
    """Metformin's Babel clique is over the cap, and CAS:1115-70-4 sorts first. As hub it anchored a
    side node of branded products; the clique's canonical id must be the hub."""
    w = ERWeights(clique_cap=3)
    ids = ["CAS:1115-70-4", "CHEBI:6801", "RXCUI:1", "RXCUI:2", "UMLS:C1"]
    ev = list(clique_evidence(ids, "babel", w, head="CHEBI:6801"))
    assert len(ev) == len(ids) - 1
    assert all("CHEBI:6801" in (a, b) for a, b, _g, _wt in ev)


def test_star_hub_falls_back_to_lexically_smallest_without_a_usable_head():
    w = ERWeights(clique_cap=3)
    ids = ["B:1", "A:1", "C:1", "D:1"]
    for head in (None, "NOT:IN_SET"):
        assert all("A:1" in (a, b) for a, b, _g, _wt in clique_evidence(ids, "babel", w, head=head))
