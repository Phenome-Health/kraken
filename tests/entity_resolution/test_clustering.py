"""Tests for connected components and Leiden/CPM clustering."""

import pytest

from kraken.entity_resolution.clustering import (
    cluster_pairs,
    connected_components,
    leiden_cpm,
)

leidenalg = pytest.importorskip("leidenalg")


def test_connected_components_deterministic():
    pairs = [("B:1", "C:1", 1.0), ("A:1", "A:2", 1.0), ("C:1", "D:1", 1.0)]
    comps = connected_components(pairs)
    assert comps == [["A:1", "A:2"], ["B:1", "C:1", "D:1"]]  # sorted, ordered by min member


def test_isolated_pair_splits_below_gamma():
    # single edge weight 0.3, gamma 0.5 -> split iff gamma > w -> split
    clusters = cluster_pairs([("A:1", "B:1", 0.3)], gamma=0.5)
    assert sorted(clusters) == [["A:1"], ["B:1"]]


def test_isolated_pair_merges_above_gamma():
    clusters = cluster_pairs([("A:1", "B:1", 0.9)], gamma=0.5)
    assert clusters == [["A:1", "B:1"]]


def test_two_cliques_weak_bridge_separate():
    # Two tight triangles joined by one weak edge; CPM should keep them apart.
    strong = 1.0
    pairs = [
        ("A:1", "A:2", strong),
        ("A:1", "A:3", strong),
        ("A:2", "A:3", strong),
        ("B:1", "B:2", strong),
        ("B:1", "B:3", strong),
        ("B:2", "B:3", strong),
        ("A:3", "B:1", 0.1),  # weak bridge
    ]
    clusters = cluster_pairs(pairs, gamma=0.5)
    clusters_sets = [set(c) for c in clusters]
    assert {"A:1", "A:2", "A:3"} in clusters_sets
    assert {"B:1", "B:2", "B:3"} in clusters_sets


def test_leiden_is_deterministic():
    pairs = [
        ("A:1", "A:2", 1.0),
        ("A:1", "A:3", 1.0),
        ("A:2", "A:3", 1.0),
        ("B:1", "B:2", 1.0),
        ("A:3", "B:1", 0.1),
    ]
    nodes = ["A:1", "A:2", "A:3", "B:1", "B:2"]
    r1 = leiden_cpm(nodes, pairs, gamma=0.5, seed=42)
    r2 = leiden_cpm(nodes, pairs, gamma=0.5, seed=42)
    assert r1 == r2


def test_higher_gamma_splits_more():
    pairs = [
        ("A:1", "A:2", 1.0),
        ("A:2", "A:3", 1.0),
        ("A:1", "A:3", 1.0),
    ]
    low = cluster_pairs(pairs, gamma=0.5)
    high = cluster_pairs(pairs, gamma=2.0)  # gamma > weight -> pull apart
    assert len(low) < len(high)
