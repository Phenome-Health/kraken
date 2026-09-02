"""Tests for name normalization and name-similarity grouping."""

from kraken.entity_resolution.name_sim import (
    group_by_normalized_name,
    is_droppable,
    name_similarity_edges,
    normalize_name,
)


def test_normalize_basic():
    assert normalize_name("  Adams-Oliver Syndrome 1 ") == "adams oliver syndrome 1"
    assert normalize_name("Café-Résumé") == "cafe resume"  # combining accents stripped
    assert normalize_name("β-amyloid") == "β amyloid"  # non-combining Greek letter kept as-is
    assert normalize_name(None) == ""
    assert normalize_name("HbA1c!!") == "hba1c"


def test_is_droppable():
    assert is_droppable("")
    assert is_droppable("ab", min_length=3)
    assert is_droppable("123")  # purely numeric
    assert is_droppable("point in time")  # stoplist
    assert not is_droppable("insulin")


def test_group_primary_names_only():
    pairs = [
        ("A:1", "Insulin"),
        ("B:1", "insulin"),
        ("C:1", "INSULIN"),
        ("D:1", "glucose"),
        ("E:1", "12345"),  # numeric -> dropped
        ("F:1", "ab"),  # too short -> dropped
    ]
    groups = group_by_normalized_name(pairs)
    assert groups == {"insulin": ["A:1", "B:1", "C:1"]}  # glucose singleton dropped


def test_name_similarity_edges_and_cap():
    groups = {"insulin": ["A:1", "B:1", "C:1"]}
    edges = list(name_similarity_edges(groups))
    assert set(edges) == {("A:1", "B:1"), ("A:1", "C:1"), ("B:1", "C:1")}
    # oversized group skipped
    big = {"x": [f"N:{i}" for i in range(50)]}
    assert list(name_similarity_edges(big, group_cap=40)) == []
