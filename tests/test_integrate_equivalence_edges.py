"""Retaining cross-cluster equivalence as close_match edges + dropping self-loops."""

import json
from types import SimpleNamespace

import jsonlines

from kraken.integrate import _EDGE_SORT_SEP, _write_equivalence_edges


def _edges_from_keyed(path):
    return [json.loads(line.split(_EDGE_SORT_SEP, 1)[1]) for line in path.read_text().splitlines()]


def test_cross_cluster_equivalence_becomes_close_match_self_loops_dropped(tmp_path):
    # Source asserts A:1 ~ {B:2, C:3}. Clustering merged A:1+B:2 (REP1) but split C:3 (REP2).
    nodes = [
        {
            "id": "A:1",
            "categories": ["biolink:Gene"],
            "provided_by": ["infores:some-source"],
            "equivalent_ids": ["A:1", "B:2", "C:3"],
        }
    ]
    src = tmp_path / "h" / "src"
    src.mkdir(parents=True)
    nodes_file = src / "nodes.jsonl"
    with jsonlines.open(nodes_file, "w") as w:
        w.write_all(nodes)

    node_map = {"A:1": "REP1", "B:2": "REP1", "C:3": "REP2"}
    config = SimpleNamespace(
        sources_to_use=["src"],
        all_harmonized_paths_resolved={"src": (nodes_file, tmp_path / "none.jsonl")},
    )
    out = tmp_path / "keyed.tsv"
    with open(out, "w") as f:
        _write_equivalence_edges(node_map, config, f)

    edges = _edges_from_keyed(out)
    # A:1~B:2 collapsed to a REP1 self-loop -> dropped; A:1~C:3 -> one REP1~REP2 close_match
    assert all(e["subject"] != e["object"] for e in edges)  # no self-loops
    cross = {(e["subject"], e["object"]) for e in edges if e["predicate"] == "biolink:close_match"}
    assert cross == {("REP1", "REP2")}  # symmetric, sorted
    edge = edges[0]
    assert edge["primary_knowledge_source"] == "infores:some-source"
    assert edge["knowledge_level"] == "knowledge_assertion"
    assert edge["agent_type"] == "not_provided"  # synthesized edge -> agent unknown
    # KRAKEN is the aggregator that produced the edge, as on every directly ingested edge
    assert edge["aggregator_knowledge_source"] == ["kraken"]


def _write_kg2(tmp_path, edges):
    src = tmp_path / "h" / "kg2"
    src.mkdir(parents=True)
    (src / "nodes.jsonl").write_text("")
    edges_file = src / "edges.jsonl"
    with jsonlines.open(edges_file, "w") as w:
        w.write_all(edges)
    return SimpleNamespace(
        sources_to_use=["kg2"],
        all_harmonized_paths_resolved={"kg2": (src / "nodes.jsonl", edges_file)},
    )


def _kg2_edge(subject, obj, pre_ids):
    return {
        "subject": subject,
        "predicate": "biolink:subclass_of",
        "object": obj,
        "primary_knowledge_source": "infores:umls-metathesaurus",
        "knowledge_level": "knowledge_assertion",
        "agent_type": "manual_agent",
        "attributes": {"infores:rtx-kg2": {"kg2pre_ids": pre_ids}},
    }


def _remapped(tmp_path, config, node_map):
    from kraken.integrate import _write_keyed_edges

    out = tmp_path / "keyed.tsv"
    _write_keyed_edges(node_map, config, out)
    return [
        (e["subject"], e["predicate"], e["object"])
        for e in _edges_from_keyed(out)
        if e["predicate"] != "biolink:close_match"
    ]


def test_kg2_originals_listed_backwards_are_re_oriented_by_cluster(tmp_path):
    """KG2 stores CHEBI:5781 subclass_of MESH:D015065 but carries both the forward assertion (RN: C0020268 ->
    C0000163) and its inverse (RB: C0000163 -> C0020268). The originals' clusters say which end is which."""
    edge = _kg2_edge(
        "CHEBI:5781",
        "MESH:D015065",
        [
            "UMLS:C0020268---UMLS:RN---None---None---None---UMLS:C0000163---src",  # listed in the stored direction
            "UMLS:C0000163---UMLS:RB---None---None---None---UMLS:C0020268---src",  # listed BACKWARDS
        ],
    )
    node_map = {
        "CHEBI:5781": "CHILD",
        "UMLS:C0020268": "CHILD",
        "MESH:D015065": "PARENT",
        "UMLS:C0000163": "PARENT",
    }
    triples = _remapped(tmp_path, _write_kg2(tmp_path, [edge]), node_map)
    assert triples == [("CHILD", "biolink:subclass_of", "PARENT")] * 2  # both come out child -> parent


def test_undetermined_orientation_follows_how_the_same_relation_oriented(tmp_path):
    """When neither original shares a cluster with a stored endpoint, the relation decides: RB pairs were seen
    to run backwards elsewhere, so this one is reversed too."""
    decided = _kg2_edge("CHEBI:1", "MESH:1", ["UMLS:P---UMLS:RB---None---None---None---UMLS:C---src"])
    undetermined = _kg2_edge("CHEBI:2", "MESH:2", ["UMLS:P2---UMLS:RB---None---None---None---UMLS:C2---src"])
    node_map = {
        "CHEBI:1": "C",
        "UMLS:C": "C",
        "MESH:1": "P",
        "UMLS:P": "P",  # decided: RB runs backwards
        "CHEBI:2": "X",
        "MESH:2": "Y",
        "UMLS:C2": "C2",
        "UMLS:P2": "P2",  # originals in clusters of their own
    }
    triples = _remapped(tmp_path, _write_kg2(tmp_path, [decided, undetermined]), node_map)
    assert ("C2", "biolink:subclass_of", "P2") in triples
    assert ("P2", "biolink:subclass_of", "C2") not in triples


def test_non_kg2_originals_are_never_re_oriented(tmp_path):
    """ROBOKOP's original_subject/original_object match its own fields, so even when an original is clustered
    with the other end, the edge keeps the recorded direction."""
    src = tmp_path / "h" / "robokop"
    src.mkdir(parents=True)
    (src / "nodes.jsonl").write_text("")
    edge = {
        "subject": "CAID:CA1",
        "predicate": "biolink:affects",
        "object": "NCBIGene:5",
        "primary_knowledge_source": "infores:gtex",
        "knowledge_level": "knowledge_assertion",
        "agent_type": "manual_agent",
        "attributes": {"infores:robokop-kg": {"original_subject": "HGVS:x", "original_object": "ENSEMBL:y"}},
    }
    with jsonlines.open(src / "edges.jsonl", "w") as w:
        w.write_all([edge])
    config = SimpleNamespace(
        sources_to_use=["robokop"],
        all_harmonized_paths_resolved={"robokop": (src / "nodes.jsonl", src / "edges.jsonl")},
    )
    # the originals landed crosswise (as a conflation elsewhere could make them)
    node_map = {"CAID:CA1": "V", "NCBIGene:5": "G", "HGVS:x": "G", "ENSEMBL:y": "V"}
    assert _remapped(tmp_path, config, node_map) == [("G", "biolink:affects", "V")]  # recorded order, not flipped


def test_babel_same_as_between_clusters_becomes_close_match(tmp_path):
    from kraken.integrate import _write_keyed_edges

    src = tmp_path / "h" / "babel"
    src.mkdir(parents=True)
    (src / "nodes.jsonl").write_text("")
    edges_file = src / "edges.jsonl"

    def edge(subject, predicate, object_):
        return {
            "subject": subject,
            "predicate": predicate,
            "object": object_,
            "primary_knowledge_source": "infores:sri-node-normalizer",
            "knowledge_level": "knowledge_assertion",
            "agent_type": "automated_agent",
        }

    with jsonlines.open(edges_file, "w") as w:
        w.write_all(
            [
                edge("MONDO:1", "biolink:same_as", "DOID:2"),  # same cluster -> self-loop
                edge("MONDO:1", "biolink:same_as", "UMLS:C3"),  # entity resolution split it off
                edge("RXCUI:1", "biolink:has_active_ingredient", "CHEBI:1"),  # not same_as -> stays as is
            ]
        )
    config = SimpleNamespace(
        sources_to_use=["babel"],
        all_harmonized_paths_resolved={"babel": (src / "nodes.jsonl", edges_file)},
    )
    node_map = {"MONDO:1": "R1", "DOID:2": "R1", "UMLS:C3": "R2", "RXCUI:1": "D", "CHEBI:1": "C"}
    out = tmp_path / "keyed.tsv"
    _write_keyed_edges(node_map, config, out)

    assert {(e["subject"], e["predicate"], e["object"]) for e in _edges_from_keyed(out)} == {
        ("R1", "biolink:close_match", "R2"),
        ("D", "biolink:has_active_ingredient", "C"),
    }


def test_any_sources_equivalence_between_clusters_becomes_close_match(tmp_path):
    """Not only Babel's: a source's same_as / exact_match that entity resolution kept apart (e.g. a kg2 same_as
    Babel overrules) can't claim the two nodes are the same either."""
    from kraken.integrate import _write_keyed_edges

    src = tmp_path / "h" / "umls"
    src.mkdir(parents=True)
    (src / "nodes.jsonl").write_text("")
    edges_file = src / "edges.jsonl"
    base = {"primary_knowledge_source": "infores:umls", "knowledge_level": "knowledge_assertion"}
    with jsonlines.open(edges_file, "w") as w:
        w.write_all(
            [
                {**base, "subject": "UMLS:C1", "predicate": "biolink:exact_match", "object": "MESH:D2"},
                {**base, "subject": "UMLS:C1", "predicate": "biolink:same_as", "object": "NCIT:C3"},
                {**base, "subject": "UMLS:C1", "predicate": "biolink:exact_match", "object": "MESH:D4"},  # one node
            ]
        )
    config = SimpleNamespace(
        sources_to_use=["umls"],
        all_harmonized_paths_resolved={"umls": (src / "nodes.jsonl", edges_file)},
    )
    node_map = {"UMLS:C1": "R1", "MESH:D2": "R2", "NCIT:C3": "R3", "MESH:D4": "R1"}
    out = tmp_path / "keyed.tsv"
    _write_keyed_edges(node_map, config, out)

    assert {(e["subject"], e["predicate"], e["object"]) for e in _edges_from_keyed(out)} == {
        ("R1", "biolink:close_match", "R2"),
        ("R1", "biolink:close_match", "R3"),
    }


def test_integration_refuses_to_run_without_babel(tmp_path):
    import pytest

    from kraken.integrate import integrate_sources

    config = SimpleNamespace(
        sources_to_use={"kg2", "robokop"},
        integrated_dir=tmp_path / "integrated",
        integrated_debug_dir=tmp_path / "integrated" / "debug",
    )
    with pytest.raises(ValueError, match="babel"):
        integrate_sources(config, biolink=None)
