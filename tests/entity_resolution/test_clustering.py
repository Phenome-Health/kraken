"""Tests for connected components and label-propagation clustering."""

import pytest

from kraken.entity_resolution.clustering import cluster_pairs, connected_components, label_propagation

pytest.importorskip("igraph")


def test_connected_components_deterministic():
    pairs = [("B:1", "C:1", 1.0), ("A:1", "A:2", 1.0), ("C:1", "D:1", 1.0)]
    comps = connected_components(pairs)
    assert comps == [["A:1", "A:2"], ["B:1", "C:1", "D:1"]]  # sorted, ordered by min member


def test_connected_pair_merges():
    # label propagation merges what's connected -> a single edge merges the pair
    clusters = cluster_pairs([("A:1", "B:1", 0.9)])
    assert clusters == [["A:1", "B:1"]]


def test_separate_components_are_separate_clusters():
    clusters = cluster_pairs([("A:1", "A:2", 1.0), ("B:1", "B:2", 1.0)])
    assert [set(c) for c in clusters] == [{"A:1", "A:2"}, {"B:1", "B:2"}]


def test_label_propagation_deterministic():
    pairs = [
        ("A:1", "A:2", 1.0),
        ("A:1", "A:3", 1.0),
        ("A:2", "A:3", 1.0),
        ("B:1", "B:2", 1.0),
        ("A:3", "B:1", 0.1),
    ]
    nodes = ["A:1", "A:2", "A:3", "B:1", "B:2"]
    r1 = label_propagation(nodes, pairs, seed=42)
    r2 = label_propagation(nodes, pairs, seed=42)
    assert r1 == r2


def test_label_propagation_no_edges_all_singletons():
    result = label_propagation(["A:1", "B:1"], [], seed=1)
    assert sorted(result) == [["A:1"], ["B:1"]]


def test_tight_clique_stays_together():
    pairs = [("A:1", "A:2", 1.0), ("A:2", "A:3", 1.0), ("A:1", "A:3", 1.0)]
    clusters = cluster_pairs(pairs)
    assert [set(c) for c in clusters] == [{"A:1", "A:2", "A:3"}]
