"""Tests for the pairwise ER scorer."""

from pathlib import Path

import pytest

from kraken.entity_resolution.eval.scorer import (
    GoldCase,
    load_gold,
    membership_from_clusters,
    score,
    score_files,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
SEED = REPO_ROOT / "config" / "entity_resolution" / "ground_truth_seed.jsonl"


def _ace_case() -> GoldCase:
    # Two true entities that a source merged: the ACE gene/protein and the disease.
    return GoldCase(
        case="conflation:MONDO:0017609",
        kind="cross-entity conflation",
        groups=(
            frozenset({"NCBIGene:1636", "HGNC:2707", "UniProtKB:P12821"}),
            frozenset({"MONDO:0017609", "orphanet:3033", "HP:0008660"}),
        ),
    )


def test_perfect_clustering():
    gc = _ace_case()
    membership = membership_from_clusters([g for g in gc.groups])
    res = score([gc], membership)
    assert res.precision == 1.0
    assert res.recall == 1.0
    assert res.f1 == 1.0
    assert res.fp == 0 and res.fn == 0
    # 2 groups of 3 -> 3+3 must-link pairs; 3*3 cannot-link pairs
    assert res.must_link_total == 6
    assert res.cannot_link_total == 9


def test_over_merge_tanks_precision():
    gc = _ace_case()
    # Everything in one cluster: all must-links satisfied, all cannot-links violated.
    all_ids = set().union(*gc.groups)
    membership = membership_from_clusters([all_ids])
    res = score([gc], membership)
    assert res.recall == 1.0  # no must-link split
    assert res.fp == 9  # every cross pair now wrongly together
    assert res.precision == pytest.approx(6 / (6 + 9))


def test_over_split_tanks_recall():
    gc = _ace_case()
    # Each id its own cluster: no cannot-link violated, but every must-link split.
    membership = membership_from_clusters([[i] for g in gc.groups for i in g])
    res = score([gc], membership)
    assert res.fn == 6
    assert res.recall == 0.0
    assert res.fp == 0
    assert res.precision == 1.0  # convention: no positive calls -> precision 1


def test_uncovered_ids_excluded_from_metrics():
    gc = _ace_case()
    # Predict only the disease group; gene ids are absent from membership.
    membership = membership_from_clusters([{"MONDO:0017609", "orphanet:3033", "HP:0008660"}])
    res = score([gc], membership)
    # gene-gene must-links (3) uncovered; disease-disease (3) covered & correct
    assert res.must_link_uncovered == 3
    assert res.tp == 3
    # all cannot-link pairs involve a missing gene id -> uncovered
    assert res.cannot_link_uncovered == 9
    assert res.fp == 0


def test_numbered_family_all_singletons_cannot_link():
    gc = GoldCase(
        case="numbered:spinocerebellar-ataxia",
        kind="numbered name family",
        groups=(frozenset({"MONDO:a"}), frozenset({"MONDO:b"}), frozenset({"MONDO:c"})),
    )
    assert gc.groups and all(len(g) == 1 for g in gc.groups)
    # correct: all apart
    apart = membership_from_clusters([["MONDO:a"], ["MONDO:b"], ["MONDO:c"]])
    res = score([gc], apart)
    assert res.cannot_link_total == 3
    assert res.tn == 3 and res.fp == 0
    assert res.must_link_total == 0
    # wrong: two merged
    merged = membership_from_clusters([["MONDO:a", "MONDO:b"], ["MONDO:c"]])
    res2 = score([gc], merged)
    assert res2.fp == 1


def test_by_kind_breakdown():
    conflation = _ace_case()
    numbered = GoldCase(
        case="numbered:x",
        kind="numbered name family (all mutually distinct diseases)",
        groups=(frozenset({"A:1"}), frozenset({"A:2"})),
    )
    membership = membership_from_clusters([g for g in conflation.groups] + [["A:1"], ["A:2"]])
    res = score([conflation, numbered], membership)
    assert "conflation" in res.by_kind
    assert "numbered_family" in res.by_kind
    assert res.by_kind["numbered_family"]["tn"] == 1


def test_contradiction_prefers_must_link(caplog):
    # Same pair asserted must_link in one case and cannot_link in another.
    c1 = GoldCase("c1", "conflation", (frozenset({"X:1", "X:2"}),))
    c2 = GoldCase("c2", "conflation", (frozenset({"X:1"}), frozenset({"X:2"})))
    membership = membership_from_clusters([{"X:1", "X:2"}])
    res = score([c1, c2], membership)
    # must_link wins -> the pair counts as a satisfied must-link, not an FP
    assert res.tp == 1
    assert res.fp == 0
    assert res.cannot_link_total == 0


def test_load_real_seed():
    assert SEED.exists(), f"seed not found at {SEED}"
    gold = load_gold(SEED)
    assert len(gold) == 28  # 14 conflation + 14 numbered
    kinds = {g.case.split(":")[0] for g in gold}
    assert "conflation" in kinds
    assert "numbered" in kinds
    # The issue #7 case is present with two groups.
    ace = next(g for g in gold if g.case == "conflation:MONDO:0017609")
    assert len(ace.groups) == 2


def test_score_files_end_to_end(tmp_path):
    gold = load_gold(SEED)
    # Build a "perfect on covered ids" membership straight from the gold groups.
    clusters = []
    for gc in gold:
        clusters.extend(list(gc.groups))
    membership = membership_from_clusters(clusters)
    mem_path = tmp_path / "membership.json"
    import json

    mem_path.write_text(json.dumps(membership))
    res = score_files(SEED, mem_path)
    assert res.precision == 1.0
    assert res.recall == 1.0
