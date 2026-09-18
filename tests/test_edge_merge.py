"""Merging a group of same-key edges into one (integrate.merge_edges).

The merge runs in a single pass. It used to fold edges in one at a time, rebuilding every list-valued property on
each fold -- quadratic in the group's size, which stalled integration for hours on one 417,750-edge group.
"""

import time

from kraken.integrate import merge_edges


def _edge(**overrides):
    edge = {
        "subject": "A:1",
        "predicate": "biolink:related_to",
        "object": "B:1",
        "primary_knowledge_source": "infores:x",
        "knowledge_level": "knowledge_assertion",
        "agent_type": "manual_agent",
    }
    edge.update(overrides)
    return edge


def test_single_edge_is_returned_as_is():
    edge = _edge(publications=["PMID:1"])
    assert merge_edges([edge]) is edge


def test_list_properties_union_in_first_seen_order():
    merged = merge_edges(
        [
            _edge(publications=["PMID:2", "PMID:1"], aggregator_knowledge_source=["infores:a"]),
            _edge(publications=["PMID:1", "PMID:3"], aggregator_knowledge_source=["infores:b", "kraken"]),
        ]
    )
    assert merged["publications"] == ["PMID:2", "PMID:1", "PMID:3"]
    assert merged["aggregator_knowledge_source"] == ["infores:a", "infores:b", "kraken"]


def test_knowledge_level_and_agent_type_take_the_first_provided_value():
    merged = merge_edges(
        [
            _edge(knowledge_level="not_provided", agent_type="manual_agent"),
            _edge(knowledge_level="prediction", agent_type="automated_agent"),
            _edge(knowledge_level="knowledge_assertion"),
        ]
    )
    assert merged["knowledge_level"] == "prediction"
    assert merged["agent_type"] == "manual_agent"


def test_attributes_merge_per_source_slot_and_union_within_it():
    merged = merge_edges(
        [
            _edge(attributes={"infores:rtx-kg2": {"kg2pre_ids": ["x"], "score": 1}}),
            _edge(attributes={"infores:rtx-kg2": {"kg2pre_ids": ["y"], "score": 2}, "infores:other": {"k": "v"}}),
        ]
    )
    assert merged["attributes"]["infores:rtx-kg2"]["kg2pre_ids"] == ["x", "y"]
    assert merged["attributes"]["infores:rtx-kg2"]["score"] == [1, 2]  # flat values inside a slot combine
    assert merged["attributes"]["infores:other"] == {"k": "v"}  # a slot only one edge has is kept as-is


def test_key_properties_and_flat_values_keep_the_first_edges():
    merged = merge_edges([_edge(publications_info=None, extra="first"), _edge(extra="second")])
    assert merged["subject"] == "A:1"
    assert merged["extra"] == "first"


def test_missing_values_never_replace_present_ones():
    merged = merge_edges([_edge(publications=["PMID:1"]), _edge()])
    assert merged["publications"] == ["PMID:1"]


def test_merge_is_linear_in_group_size():
    """Doubling a group of all-distinct values must roughly double the work, not quadruple it."""

    def group(n):
        return [_edge(attributes={"infores:rtx-kg2": {"kg2pre_ids": [f"P:{i}"]}}) for i in range(n)]

    def timed(n):
        g = group(n)
        start = time.perf_counter()
        merge_edges(g)
        return time.perf_counter() - start

    small, large = min(timed(20_000) for _ in range(3)), min(timed(80_000) for _ in range(3))
    assert large < small * 8, f"4x the edges took {large / small:.1f}x as long -- the merge has gone superlinear"
