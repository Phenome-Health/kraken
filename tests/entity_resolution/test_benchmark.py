"""The ER benchmark runner: case outcomes from a build's membership."""

import json

import pytest

from kraken.entity_resolution.eval.benchmark import (
    DEFAULT_BENCHMARK_PATH,
    BenchmarkCase,
    evaluate,
    load_cases,
    membership_from_nodes_file,
)

T2D = BenchmarkCase(
    case="t2d",
    kind="disease",
    clusters={"type 2 diabetes": ["MONDO:0005148", "OMIM:MTHU069260"], "type 1 diabetes": ["MONDO:0005147"]},
)


def test_a_case_passes_when_clusters_are_whole_and_apart():
    result = evaluate(T2D, {"MONDO:0005148": "N1", "OMIM:MTHU069260": "N1", "MONDO:0005147": "N2"})
    assert result.passed and result.outcome == "pass"


def test_a_split_a_conflation_and_a_missing_id_each_fail_a_case():
    split = evaluate(T2D, {"MONDO:0005148": "N1", "OMIM:MTHU069260": "N3", "MONDO:0005147": "N2"})
    assert split.splits == {"type 2 diabetes": {"N1": ["MONDO:0005148"], "N3": ["OMIM:MTHU069260"]}}
    assert split.outcome == "REGRESSION"

    conflated = evaluate(T2D, {"MONDO:0005148": "N1", "OMIM:MTHU069260": "N1", "MONDO:0005147": "N1"})
    assert conflated.conflations == {
        "N1": {"type 2 diabetes": ["MONDO:0005148", "OMIM:MTHU069260"], "type 1 diabetes": ["MONDO:0005147"]}
    }

    missing = evaluate(T2D, {"MONDO:0005148": "N1", "MONDO:0005147": "N2"})
    assert missing.missing == ["OMIM:MTHU069260"] and not missing.passed


def test_a_known_issue_is_fixed_once_it_passes():
    known = BenchmarkCase(**{**T2D.__dict__, "known_issue": "OMIM:MTHU069260 is a separate node"})
    assert evaluate(known, {"MONDO:0005148": "N1", "OMIM:MTHU069260": "N3", "MONDO:0005147": "N2"}).outcome == (
        "known issue"
    )
    assert evaluate(known, {"MONDO:0005148": "N1", "OMIM:MTHU069260": "N1", "MONDO:0005147": "N2"}).outcome == "fixed"


def test_membership_comes_from_each_nodes_equivalent_ids(tmp_path):
    nodes = tmp_path / "nodes.jsonl"
    nodes.write_text(
        json.dumps({"id": "MONDO:0005148", "equivalent_ids": ["MONDO:0005148", "OMIM:MTHU069260"]})
        + "\n"
        + json.dumps({"id": "MONDO:0005147"})
        + "\n"
    )
    wanted = {"OMIM:MTHU069260", "MONDO:0005147", "NOT:THERE"}
    assert membership_from_nodes_file(nodes, wanted) == {
        "OMIM:MTHU069260": "MONDO:0005148",
        "MONDO:0005147": "MONDO:0005147",
    }


@pytest.mark.skipif(not DEFAULT_BENCHMARK_PATH.exists(), reason="no benchmark case file")
def test_the_benchmark_case_file_is_well_formed():
    cases = load_cases(DEFAULT_BENCHMARK_PATH)
    assert cases
    for case in cases:
        ids = [curie for members in case.clusters.values() for curie in members]
        assert len(ids) == len(set(ids)), f"{case.case}: an id is in two clusters"
        assert case.flagged_ids <= set(ids), f"{case.case}: a flagged id is in no cluster"


@pytest.mark.skipif(not DEFAULT_BENCHMARK_PATH.exists(), reason="no benchmark case file")
def test_no_two_cases_disagree_about_a_pair():
    """Two ids one case keeps together, no other case may keep apart."""
    from itertools import combinations

    together, apart = set(), set()
    for case in load_cases(DEFAULT_BENCHMARK_PATH):
        groups = [sorted(members) for members in case.clusters.values()]
        for members in groups:
            together.update(combinations(members, 2))
        for first, second in combinations(groups, 2):
            apart.update((a, b) if a < b else (b, a) for a in first for b in second)
    assert not together & apart, f"pairs both together and apart: {sorted(together & apart)[:5]}"


def test_comparing_two_builds_finds_what_a_case_outcome_cannot():
    """A case carrying a known_issue reads "known issue" however it fails, so breakage arriving inside an
    already-failing case never shows up as a REGRESSION. Comparing pair by pair against a baseline does."""
    from kraken.entity_resolution.eval.benchmark import compare

    case = BenchmarkCase(
        case="two-things",
        kind="chemical",
        clusters={"a thing": ["X:1", "X:2"], "another thing": ["Y:1"]},
        known_issue="X:2 is in a node of its own",
    )
    before = {"X:1": "n1", "X:2": "n1", "Y:1": "n2"}
    after = {"X:1": "n1", "X:2": "n3", "Y:1": "n1"}  # split X, and pulled Y in with X:1
    changes = compare([case], before, after)
    assert changes["two-things"]["lost"] == [("X:1", "X:2")]
    assert changes["two-things"]["gained"] == [("X:1", "Y:1")]
    assert compare([case], before, before) == {}, "nothing changed, nothing reported"
