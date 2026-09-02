"""Test the out-of-core entity-resolution build (evidence -> clusters -> nodes)."""

from types import SimpleNamespace

import jsonlines
import pytest

from kraken.entity_resolution.build import _original_endpoints, resolve_entities

pytest.importorskip("leidenalg")


def test_original_endpoints_uncanonicalizes_kg2():
    # KG2 stores canonicalized subject/object; the originals are in kg2pre_ids
    # (subject at field 0, object at field 5), and a merged edge can carry several.
    kg2_edge = {
        "subject": "UNII:1",
        "object": "PUBCHEM.COMPOUND:1",
        "predicate": "biolink:close_match",
        "attributes": {
            "infores:rtx-kg2": {
                "kg2pre_ids": [
                    "ATC:X---UMLS:xref---None---None---None---UMLS:Y---umls_source:ATC",
                    "ATC:Z---UMLS:xref---None---None---None---UMLS:Y---umls_source:ATC",
                ]
            }
        },
    }
    assert _original_endpoints(kg2_edge, "kg2") == [("ATC:X", "UMLS:Y"), ("ATC:Z", "UMLS:Y")]
    # a native (non-canonicalized) source uses its own endpoints
    assert _original_endpoints({"subject": "A:1", "object": "B:1"}, "refmet") == [("A:1", "B:1")]
    # a canonicalized aggregator without a known un-canonicalizer is skipped (not junk)
    assert _original_endpoints({"subject": "A:1", "object": "B:1"}, "robokop") is None


def test_kg2_close_match_downweighted_by_subclass_count(tmp_path):
    import io

    from kraken.entity_resolution.build import _subclass_penalized_weight, _write_kg2_match_evidence
    from kraken.entity_resolution.weights import ERWeights

    # pure penalty: each hierarchical edge multiplies by the decay
    assert _subclass_penalized_weight(0.2, 0, 0.5) == 0.2
    assert _subclass_penalized_weight(0.2, 1, 0.5) == 0.1
    assert _subclass_penalized_weight(0.2, 2, 0.5) == 0.05

    w = ERWeights()

    def edge(pred, subj, obj):
        return {
            "subject": "UNII:1",
            "object": "PUBCHEM.COMPOUND:1",
            "predicate": pred,
            "attributes": {"infores:rtx-kg2": {"kg2pre_ids": [f"{subj}---rel---None---None---None---{obj}---src"]}},
        }

    # one close_match + two subclass edges (either direction) between ATC:X and UMLS:Y
    edges = [
        edge("biolink:close_match", "ATC:X", "UMLS:Y"),
        edge("biolink:subclass_of", "ATC:X", "UMLS:Y"),
        edge("biolink:subclass_of", "UMLS:Y", "ATC:X"),
    ]
    ef = tmp_path / "kg2_edges.jsonl"
    with jsonlines.open(ef, "w") as wr:
        wr.write_all(edges)

    out = io.StringIO()
    _write_kg2_match_evidence(ef, w, out)
    lines = [ln for ln in out.getvalue().splitlines() if ln]
    assert len(lines) == 1
    a, b, _group, weight = lines[0].split("\t")
    assert (a, b) == ("ATC:X", "UMLS:Y")
    assert float(weight) == w.close_match_weight * (w.subclass_penalty_decay**2)  # 2 hierarchical edges


def _write_source(tmp_path, source, records):
    d = tmp_path / "harmonized" / source
    d.mkdir(parents=True)
    nodes_path = d / "nodes.jsonl"
    edges_path = d / "edges.jsonl"
    with jsonlines.open(nodes_path, "w") as w:
        w.write_all(records)
    return nodes_path, edges_path


def _fake_config(tmp_path):
    gene_clique = ["NCBIGene:1636", "HGNC:2707", "UniProtKB:P12821"]
    harmonized = {
        "ncbigene": _write_source(
            tmp_path,
            "ncbigene",
            [
                {
                    "id": "NCBIGene:1636",
                    "categories": ["biolink:Gene"],
                    "provided_by": ["ncbigene"],
                    "equivalent_ids": gene_clique,
                    "name": "ACE",
                },
                {
                    "id": "HGNC:2707",
                    "categories": ["biolink:Gene"],
                    "provided_by": ["ncbigene"],
                    "equivalent_ids": gene_clique,
                },
                {
                    "id": "UniProtKB:P12821",
                    "categories": ["biolink:Protein"],
                    "provided_by": ["ncbigene"],
                    "equivalent_ids": gene_clique,
                },
            ],
        ),
        "umls": _write_source(
            tmp_path,
            "umls",
            [
                {
                    "id": "MONDO:0017609",
                    "categories": ["biolink:Disease"],
                    "provided_by": ["umls"],
                    "equivalent_ids": ["MONDO:0017609", "orphanet:3033"],
                    "name": "renal tubular dysgenesis",
                },
                {
                    "id": "orphanet:3033",
                    "categories": ["biolink:Disease"],
                    "provided_by": ["umls"],
                    "equivalent_ids": ["MONDO:0017609", "orphanet:3033"],
                },
            ],
        ),
        # weak aggregator conflation across gene and disease
        "kg2": _write_source(
            tmp_path,
            "kg2",
            [
                {
                    "id": "NCBIGene:1636",
                    "categories": [],
                    "provided_by": ["kg2"],
                    "equivalent_ids": ["NCBIGene:1636", "MONDO:0017609"],
                },
            ],
        ),
    }
    integrated = tmp_path / "integrated"
    return SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )


def test_build_splits_conflation_and_writes_nodes(tmp_path):
    config = _fake_config(tmp_path)
    node_id_to_rep = resolve_entities(config, biolink=None)

    # gene/protein share a representative; disease is separate
    assert node_id_to_rep["NCBIGene:1636"] == node_id_to_rep["HGNC:2707"] == node_id_to_rep["UniProtKB:P12821"]
    assert node_id_to_rep["MONDO:0017609"] == node_id_to_rep["orphanet:3033"]
    assert node_id_to_rep["NCBIGene:1636"] != node_id_to_rep["MONDO:0017609"]

    # representative id for the gene cluster is HGNC (prefix ranking)
    assert node_id_to_rep["NCBIGene:1636"] == "HGNC:2707"

    # canonical nodes file: one node per cluster, gene node carries all three ids
    nodes = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes.extend(r)
    gene_node = next(n for n in nodes if n["id"] == "HGNC:2707")
    assert set(gene_node["equivalent_ids"]) >= {"NCBIGene:1636", "HGNC:2707", "UniProtKB:P12821"}
    assert gene_node["name"] == "ACE"
    disease_node = next(n for n in nodes if n["id"] == "MONDO:0017609")
    assert "orphanet:3033" in disease_node["equivalent_ids"]


def test_canonical_equiv_ids_are_disjoint_and_obey_one_id(tmp_path):
    # A source node that lists TWO LM ids as equivalent conflates two structural
    # entities. The one-LM-per-cluster guardrail must split them, and no canonical
    # node may carry both (equivalent_ids come from cluster membership, not the
    # source's raw list).
    harmonized = {
        "lipidmaps": _write_source(
            tmp_path,
            "lipidmaps",
            [
                {
                    "id": "LM:1",
                    "categories": ["biolink:SmallMolecule"],
                    "provided_by": ["lipidmaps"],
                    "equivalent_ids": ["LM:1", "LM:2"],
                    "name": "lipid one",
                },
                {
                    "id": "LM:2",
                    "categories": ["biolink:SmallMolecule"],
                    "provided_by": ["lipidmaps"],
                    "equivalent_ids": ["LM:1", "LM:2"],
                    "name": "lipid two",
                },
            ],
        ),
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    resolve_entities(config, biolink=None)
    nodes = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes.extend(r)
    # no canonical node has two LM ids
    for n in nodes:
        assert sum(e.startswith("LM:") for e in n["equivalent_ids"]) <= 1
    # equivalent_ids are disjoint across canonical nodes
    seen: set[str] = set()
    for n in nodes:
        eqs = set(n["equivalent_ids"])
        assert not (eqs & seen), f"overlapping equivalent_ids: {eqs & seen}"
        seen |= eqs


def test_bare_id_nn_label_is_retained_as_synonym(tmp_path, monkeypatch):
    # A bare equivalency-list id (no harmonized node) gets its name from the node
    # normalizer; that name must survive as a synonym on the merged node.
    from kraken.entity_resolution import build as build_mod
    from kraken.entity_resolution.sri_nodenorm import NormInfo

    class _StubNN:
        def __init__(self, *a, **k):
            pass

        def resolve(self, curies, **k):
            return {
                c: NormInfo(label="Angiotensin Converting Enzyme", categories=("biolink:Gene",))
                for c in curies
                if c == "MESH:D000806"
            }

        def close(self):
            pass

    monkeypatch.setattr(build_mod, "NodeNormClient", _StubNN)

    clique = ["NCBIGene:1636", "HGNC:2707", "MESH:D000806"]  # MESH id is bare
    harmonized = {
        "ncbigene": _write_source(
            tmp_path,
            "ncbigene",
            [
                {
                    "id": "NCBIGene:1636",
                    "categories": ["biolink:Gene"],
                    "provided_by": ["ncbigene"],
                    "equivalent_ids": clique,
                    "name": "ACE",
                },
                {
                    "id": "HGNC:2707",
                    "categories": ["biolink:Gene"],
                    "provided_by": ["ncbigene"],
                    "equivalent_ids": clique,
                },
            ],
        ),
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    resolve_entities(config, biolink=None)
    nodes = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes.extend(r)
    node = next(n for n in nodes if "NCBIGene:1636" in n["equivalent_ids"])
    assert "MESH:D000806" in node["equivalent_ids"]  # bare id joined the cluster
    assert "Angiotensin Converting Enzyme" in node.get("synonyms", [])  # NN label kept as synonym


def test_name_only_pair_merges_when_compatible(tmp_path):
    # Two nodes linked ONLY by a shared normalized name, same branch, no enforced
    # id -> they should merge (name weight is above gamma).
    harmonized = {
        "src": _write_source(
            tmp_path,
            "src",
            [
                {
                    "id": "FOO:1",
                    "categories": ["biolink:SmallMolecule"],
                    "provided_by": ["src"],
                    "equivalent_ids": ["FOO:1"],
                    "name": "Shared Compound Name",
                },
                {
                    "id": "FOO:2",
                    "categories": ["biolink:SmallMolecule"],
                    "provided_by": ["src"],
                    "equivalent_ids": ["FOO:2"],
                    "name": "shared  compound   name",
                },  # normalizes identically
            ],
        ),
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    m = resolve_entities(config, biolink=None)
    assert m["FOO:1"] == m["FOO:2"]  # merged on the shared name alone


def test_name_collision_across_branches_does_not_merge(tmp_path):
    # Same shared name, but one is a gene and one a disease -> the branch guardrail
    # prunes the name edge at formation, so they stay separate.
    harmonized = {
        "g": _write_source(
            tmp_path,
            "g",
            [
                {
                    "id": "GENE:1",
                    "categories": ["biolink:Gene"],
                    "provided_by": ["g"],
                    "equivalent_ids": ["GENE:1"],
                    "name": "insulin",
                },
            ],
        ),
        "d": _write_source(
            tmp_path,
            "d",
            [
                {
                    "id": "DIS:1",
                    "categories": ["biolink:Disease"],
                    "provided_by": ["d"],
                    "equivalent_ids": ["DIS:1"],
                    "name": "Insulin",
                },
            ],
        ),
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    m = resolve_entities(config, biolink=None)
    assert m["GENE:1"] != m["DIS:1"]  # pruned: gene vs disease never merges on a name


def test_kg2_match_edge_merges_original_endpoints_not_canonical(tmp_path):
    # A KG2 exact_match edge stored on canonical endpoints (UNII/PUBCHEM) must be
    # un-canonicalized to its originals (ATC:X = UMLS:Y) before entering the match
    # graph, so the originals merge and the canonical ids are never used as evidence.
    src = _write_source(
        tmp_path,
        "src",
        [
            {
                "id": "ATC:X",
                "categories": ["biolink:SmallMolecule"],
                "provided_by": ["src"],
                "equivalent_ids": ["ATC:X"],
                "name": "chem x",
            },
            {
                "id": "UMLS:Y",
                "categories": ["biolink:SmallMolecule"],
                "provided_by": ["src"],
                "equivalent_ids": ["UMLS:Y"],
                "name": "chem y",
            },
        ],
    )
    kg2_nodes, kg2_edges = _write_source(tmp_path, "kg2", [])  # no kg2 nodes, just the edge
    with jsonlines.open(kg2_edges, "w") as w:
        w.write(
            {
                "subject": "UNII:1",
                "object": "PUBCHEM.COMPOUND:1",
                "predicate": "biolink:exact_match",
                "primary_knowledge_source": "infores:atc-codes-umls",
                "knowledge_level": "knowledge_assertion",
                "agent_type": "manual_agent",
                "attributes": {
                    "infores:rtx-kg2": {
                        "kg2pre_ids": ["ATC:X---UMLS:xref---None---None---None---UMLS:Y---umls_source:ATC"]
                    }
                },
            }
        )
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"src": src, "kg2": (kg2_nodes, kg2_edges)},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    m = resolve_entities(config, biolink=None)
    assert m["ATC:X"] == m["UMLS:Y"]  # merged via the un-canonicalized exact_match
    assert "UNII:1" not in m and "PUBCHEM.COMPOUND:1" not in m  # canonical ids never used


def test_build_isolated_node_becomes_singleton(tmp_path):
    # a node with no equivalencies and a unique name is its own cluster
    harmonized = {
        "refmet": _write_source(
            tmp_path,
            "refmet",
            [
                {
                    "id": "RM:1",
                    "categories": ["biolink:SmallMolecule"],
                    "provided_by": ["refmet"],
                    "equivalent_ids": ["RM:1"],
                    "name": "some unique metabolite name",
                },
            ],
        ),
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    node_id_to_rep = resolve_entities(config, biolink=None)
    assert node_id_to_rep["RM:1"] == "RM:1"
    nodes = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes.extend(r)
    assert len(nodes) == 1
    assert nodes[0]["id"] == "RM:1"
