"""Retaining cross-cluster equivalence as same_as edges + dropping self-loops."""

import json
from types import SimpleNamespace

import jsonlines

from kraken.integrate import _EDGE_SORT_SEP, _write_equivalence_edges


def _edges_from_keyed(path):
    return [json.loads(line.split(_EDGE_SORT_SEP, 1)[1]) for line in path.read_text().splitlines()]


def test_cross_cluster_equivalence_becomes_same_as_self_loops_dropped(tmp_path):
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
        er_nodenorm_cache_path=tmp_path / "nn.sqlite",  # empty cache -> no NN cliques
    )
    out = tmp_path / "keyed.tsv"
    with open(out, "w") as f:
        _write_equivalence_edges(node_map, config, f)

    edges = _edges_from_keyed(out)
    # A:1~B:2 collapsed to a REP1 self-loop -> dropped; A:1~C:3 -> one REP1~REP2 same_as
    assert all(e["subject"] != e["object"] for e in edges)  # no self-loops
    same_as = {(e["subject"], e["object"]) for e in edges if e["predicate"] == "biolink:same_as"}
    assert same_as == {("REP1", "REP2")}  # symmetric, sorted
    edge = edges[0]
    assert edge["primary_knowledge_source"] == "infores:some-source"
    assert edge["knowledge_level"] == "knowledge_assertion"
    assert edge["agent_type"] == "not_provided"  # synthesized edge -> agent unknown


def test_nn_clique_equivalence_edge_uses_normalizer_primary_ks(tmp_path):
    from kraken.entity_resolution.sri_nodenorm import NodeNormClient, NormInfo

    # Seed a normalizer clique into the cache: MONDO:1 (canonical) ~ DOID:2.
    cache = tmp_path / "nn.sqlite"
    client = NodeNormClient(cache)
    info = NormInfo(label="d", categories=("biolink:Disease",), canonical="MONDO:1")
    client._cache_put("MONDO:1", info, resolved=True)
    client._cache_put("DOID:2", info, resolved=True)
    client.close()

    empty_nodes = tmp_path / "empty.jsonl"
    empty_nodes.write_text("")
    node_map = {"MONDO:1": "MONDO:1", "DOID:2": "DOID:2"}  # split into different clusters
    config = SimpleNamespace(
        sources_to_use=["src"],
        all_harmonized_paths_resolved={"src": (empty_nodes, tmp_path / "none.jsonl")},
        er_nodenorm_cache_path=cache,
    )
    out = tmp_path / "keyed.tsv"
    with open(out, "w") as f:
        _write_equivalence_edges(node_map, config, f)

    edges = _edges_from_keyed(out)
    assert len(edges) == 1
    assert (edges[0]["subject"], edges[0]["object"]) == ("DOID:2", "MONDO:1")  # sorted
    assert edges[0]["primary_knowledge_source"] == "infores:sri-node-normalizer"
