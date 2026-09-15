"""Test the out-of-core entity-resolution build (evidence -> clusters -> nodes)."""

from types import SimpleNamespace

import jsonlines
import pytest

from kraken.entity_resolution.build import _original_endpoints, resolve_entities

pytest.importorskip("igraph")


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
    a, b, _group, weight, kind = lines[0].split("\t")
    assert kind == "match:kg2"
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


def test_babel_label_of_a_listed_id_is_retained_as_synonym(tmp_path):
    # An id a source only lists (MESH:D000806) is named by its own Babel node; that name must survive as a synonym
    # on the merged node.
    clique = ["NCBIGene:1636", "HGNC:2707", "MESH:D000806"]
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
        "babel": _write_source_with_edges(
            tmp_path,
            "babel",
            [
                {
                    "id": "MESH:D000806",
                    "categories": ["biolink:Gene"],
                    "provided_by": ["infores:sri-node-normalizer"],
                    "name": "Angiotensin Converting Enzyme",
                }
            ],
            [],
        ),
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    resolve_entities(config, biolink=None)
    nodes = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes.extend(r)
    node = next(n for n in nodes if "NCBIGene:1636" in n["equivalent_ids"])
    assert "MESH:D000806" in node["equivalent_ids"]  # the listed id joined the cluster
    assert "Angiotensin Converting Enzyme" in node.get("synonyms", [])  # its Babel label kept as a synonym


def test_name_only_pair_merges_when_compatible(tmp_path):
    # Two nodes linked ONLY by a shared normalized name, same branch, no enforced
    # id -> they should merge (name weight reaches tau).
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
    )
    m = resolve_entities(config, biolink=None)
    assert m["ATC:X"] == m["UMLS:Y"]  # merged via the un-canonicalized exact_match
    assert "UNII:1" not in m and "PUBCHEM.COMPOUND:1" not in m  # canonical ids never used


def test_inherited_cats_from_single_family_referencing_nodes(tmp_path):
    # A bare id (only in an equiv list) inherits its category from the single-family
    # node that lists it; a conflated multi-family node does NOT propagate.
    from kraken.entity_resolution.build import _stage1_write_evidence_and_facts
    from kraken.entity_resolution.families import BranchFamilies
    from kraken.entity_resolution.id_facts import IdFactsStore
    from kraken.entity_resolution.weights import ERWeights

    fams = BranchFamilies.load()
    src, _edges = _write_source(
        tmp_path,
        "src",
        [
            {
                "id": "UBERON:1",
                "categories": ["biolink:AnatomicalEntity"],
                "provided_by": ["src"],
                "equivalent_ids": ["UBERON:1", "AEO:1"],
            },  # clean single-family -> propagates
            {
                "id": "KG2NODE:1",
                "categories": ["biolink:Disease", "biolink:Gene"],
                "provided_by": ["kg2"],
                "equivalent_ids": ["KG2NODE:1", "BARE:2"],
            },  # conflated multi-family -> does NOT propagate
        ],
    )
    config = SimpleNamespace(all_harmonized_paths_resolved={"src": (src, tmp_path / "none.jsonl")})
    facts = IdFactsStore(tmp_path / "facts.sqlite")
    inherited, _taxon, node_ids, _seeds, _prov = _stage1_write_evidence_and_facts(
        config, ERWeights(), fams, tmp_path / "ev.tmp", tmp_path / "nm.tmp", {"src": 0}, facts
    )
    facts.close()
    assert inherited["AEO:1"] == {"biolink:AnatomicalEntity"}  # bare id typed via its referencing node
    assert inherited["UBERON:1"] == {"biolink:AnatomicalEntity"}  # self-typed
    assert "BARE:2" not in inherited  # conflated multi-family node did not spread its typing
    assert node_ids == {"UBERON:1", "KG2NODE:1"}


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
    )
    node_id_to_rep = resolve_entities(config, biolink=None)
    assert node_id_to_rep["RM:1"] == "RM:1"
    nodes = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes.extend(r)
    assert len(nodes) == 1
    assert nodes[0]["id"] == "RM:1"


def _write_source_with_edges(tmp_path, source, nodes, edges):
    """A harmonized source with both files written (``_write_source`` writes nodes only)."""
    d = tmp_path / "harmonized" / source
    d.mkdir(parents=True)
    nodes_path, edges_path = d / "nodes.jsonl", d / "edges.jsonl"
    with jsonlines.open(nodes_path, "w") as w:
        w.write_all(nodes)
    with jsonlines.open(edges_path, "w") as w:
        w.write_all(edges)
    return nodes_path, edges_path


def test_aggregator_original_endpoints_become_resolvable_ids(tmp_path):
    """An aggregator edge's ORIGINAL endpoint resolves even when it is in no node set.

    This is robokop's gtex shape: the edge is stored on a ``CAID:`` variant node but originates at
    an ``HGVS:`` expression that exists nowhere as a node and that Babel doesn't know.
    Edges are remapped by their originals, so before originals were treated as real ids every one
    of these edges was dropped as an orphan (18.5M of them in the 2.1.1 build).
    """
    hgvs = "HGVS:NC_000021.9:g.25840043C>G"
    nodes = [
        {
            "id": "CAID:CA15984545",
            "categories": ["biolink:SequenceVariant"],
            "provided_by": ["infores:robokop-kg"],
            "equivalent_ids": ["CAID:CA15984545", "DBSNP:rs1827747"],
            "name": "rs1827747",
        },
        {
            "id": "NCBIGene:2551",
            "categories": ["biolink:Gene"],
            "provided_by": ["infores:robokop-kg"],
            "equivalent_ids": ["NCBIGene:2551"],
            "name": "GABPA",
        },
    ]
    edges = [
        {
            "subject": "CAID:CA15984545",
            "object": "NCBIGene:2551",
            "predicate": "biolink:affects",
            "primary_knowledge_source": "infores:gtex",
            "attributes": {
                "infores:robokop-kg": {"original_subject": hgvs, "original_object": "ENSEMBL:ENSG00000154727"}
            },
        }
    ]
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"robokop": _write_source_with_edges(tmp_path, "robokop", nodes, edges)},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    node_id_to_rep = resolve_entities(config, biolink=None)

    # The original ids are resolvable -- so the edge can be remapped instead of orphaned...
    assert hgvs in node_id_to_rep
    assert "ENSEMBL:ENSG00000154727" in node_id_to_rep
    # ...and they land on the node the aggregator said they were, not on duplicates of it.
    assert node_id_to_rep[hgvs] == node_id_to_rep["CAID:CA15984545"]
    assert node_id_to_rep["ENSEMBL:ENSG00000154727"] == node_id_to_rep["NCBIGene:2551"]

    # The HGVS id inherits the variant's category rather than falling through to NamedThing
    # (nothing else can type it: no node carries it and Babel doesn't know HGVS).
    nodes_out = []
    with jsonlines.open(config.integrated_nodes_path) as r:
        nodes_out.extend(r)
    variant = next(n for n in nodes_out if hgvs in n.get("equivalent_ids", []))
    assert variant["categories"] == ["biolink:SequenceVariant"]


def test_native_source_edges_contribute_no_aliases(tmp_path):
    """A native source's endpoints are already its own ids, so no alias evidence is created."""
    from kraken.entity_resolution.uncanonicalize import original_alias_pairs

    edge = {"subject": "A:1", "object": "B:2", "predicate": "biolink:affects"}
    assert original_alias_pairs(edge, "ncbigene") == []
    # ...and an aggregator edge with no recorded originals contributes none either.
    assert original_alias_pairs(edge, "robokop") == []


def test_kg2_parallel_close_matches_from_different_kses_reach_tau_through_the_real_writer(tmp_path):
    """The primary KS has to survive kg2's two-phase buffered writer, or parallel claims collapse."""
    import io

    from kraken.entity_resolution.build import _write_kg2_match_evidence
    from kraken.entity_resolution.match_graph import accumulate
    from kraken.entity_resolution.weights import ERWeights

    w = ERWeights()

    def edge(ks):
        return {
            "subject": "UNII:1",
            "object": "PUBCHEM.COMPOUND:1",
            "predicate": "biolink:close_match",
            "primary_knowledge_source": ks,
            "attributes": {"infores:rtx-kg2": {"kg2pre_ids": ["ATC:X---rel---None---None---None---UMLS:Y---src"]}},
        }

    ef = tmp_path / "kg2_edges.jsonl"
    with jsonlines.open(ef, "w") as wr:
        wr.write_all([edge("infores:mesh"), edge("infores:go"), edge("infores:chv-umls")])

    out = io.StringIO()
    _write_kg2_match_evidence(ef, w, out)
    rows = [ln.split("\t")[:4] for ln in out.getvalue().splitlines() if ln]
    assert len({group for _a, _b, group, _wt in rows}) == 3, "each KS should land in its own group"
    total = accumulate([(a, b, g, float(wt)) for a, b, g, wt in rows], w)[("ATC:X", "UMLS:Y")]
    assert total >= w.tau, f"three independent KSes should reach tau, got {total}"


def test_alias_is_skipped_when_its_original_is_already_in_the_graph(tmp_path):
    """An alias only attaches ids nothing else places. If any source's nodes or lists already carry the original
    -- here a source that sorts AFTER robokop, so this also checks the node pass completes before aliases -- the
    alias contributes nothing, and the original stays wherever that source's evidence puts it."""
    robokop_nodes = [
        {
            "id": "CAID:CA1",
            "categories": ["biolink:SequenceVariant"],
            "provided_by": ["infores:robokop-kg"],
            "equivalent_ids": ["CAID:CA1"],
        },
        {
            "id": "NCBIGene:5",
            "categories": ["biolink:Gene"],
            "provided_by": ["infores:robokop-kg"],
            "equivalent_ids": ["NCBIGene:5"],
        },
    ]
    robokop_edges = [
        {
            "subject": "CAID:CA1",
            "object": "NCBIGene:5",
            "predicate": "biolink:affects",
            "primary_knowledge_source": "infores:gtex",
            "attributes": {
                "infores:robokop-kg": {"original_subject": "HGVS:already", "original_object": "ENSEMBL:new"}
            },
        }
    ]
    umls_nodes = [
        {
            "id": "HGVS:already",
            "categories": ["biolink:SequenceVariant"],
            "provided_by": ["infores:umls"],
            "equivalent_ids": ["HGVS:already"],
        },
    ]
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={
            "robokop": _write_source_with_edges(tmp_path, "robokop", robokop_nodes, robokop_edges),
            "umls": _write_source_with_edges(tmp_path, "umls", umls_nodes, []),
        },
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    m = resolve_entities(config, biolink=None)
    assert m["HGVS:already"] != m["CAID:CA1"], "an alias merged an id another source already places"
    assert m["ENSEMBL:new"] == m["NCBIGene:5"], "an id nothing else places should still attach via its alias"


def test_oversized_cluster_evidence_report_counts_kinds_inside_the_cluster(tmp_path, monkeypatch):
    """The kind column is what tells a Babel clique from an aggregator list; only rows linking two members of the
    same oversized cluster count. The breakdown goes to a file, one JSON line per cluster, largest first."""
    import json

    from kraken.entity_resolution import build as build_mod

    monkeypatch.setattr(build_mod, "OVERSIZED_EVIDENCE_REPORT_MIN_SIZE", 3)
    evidence = tmp_path / "ev.tmp"
    evidence.write_text(
        build_mod._evidence_row("A:1", "A:2", "babel_derived", 0.5, kind="babel")
        + build_mod._evidence_row("A:2", "A:3", "babel_derived", 0.5, kind="equiv:kg2")
        + build_mod._evidence_row("A:1", "A:3", "babel_derived", 0.5, kind="equiv:kg2")
        + build_mod._evidence_row("A:1", "B:1", "name_similarity", 0.7, kind="name_sim")  # crosses clusters
    )
    curie_to_cluster = {"A:1": 0, "A:2": 0, "A:3": 0, "B:1": 1}
    report = tmp_path / "debug" / "oversized_clusters.jsonl"
    build_mod._report_oversized_cluster_evidence(evidence, curie_to_cluster, report)
    rows = [json.loads(line) for line in report.read_text().splitlines()]
    assert len(rows) == 1  # cluster 1 is too small to report
    assert rows[0]["members"] == 3
    assert rows[0]["evidence_rows_inside_by_kind"] == {"equiv:kg2": 2, "babel": 1}
    assert rows[0]["id_prefixes"] == {"A": 3}


def test_no_report_file_when_no_cluster_is_oversized(tmp_path):
    from kraken.entity_resolution import build as build_mod

    evidence = tmp_path / "ev.tmp"
    evidence.write_text(build_mod._evidence_row("A:1", "A:2", "g", 0.5, kind="babel"))
    report = tmp_path / "oversized_clusters.jsonl"
    build_mod._report_oversized_cluster_evidence(evidence, {"A:1": 0, "A:2": 0}, report)
    assert not report.exists()


def test_babel_supplies_cliques_gene_protein_conflations_and_per_id_categories(tmp_path):
    """Babel's cliques and gene/protein conflations come from its edges, each id's category from its node. Its
    drug/chemical relations are not merge evidence (that conflation is off)."""
    from kraken.entity_resolution import build as build_mod

    def node(curie, category, name=None):
        return {"id": curie, "categories": [category], "provided_by": ["infores:sri-node-normalizer"], "name": name}

    def edge(subject, predicate, object_):
        return {
            "subject": subject,
            "predicate": predicate,
            "object": object_,
            "primary_knowledge_source": "infores:sri-node-normalizer",
        }

    babel = _write_source_with_edges(
        tmp_path,
        "babel",
        [
            node("NCBIGene:1636", "biolink:Gene", "ACE"),
            node("HGNC:2707", "biolink:Gene", "ACE gene symbol"),
            node("UniProtKB:P12821", "biolink:Protein", "Angiotensin-converting enzyme"),
            node("CHEBI:6801", "biolink:SmallMolecule", "metformin"),
            node("RXCUI:861007", "biolink:Drug", "metformin 500 MG Oral Tablet"),
        ],
        [
            edge("NCBIGene:1636", "biolink:same_as", "HGNC:2707"),
            {
                **edge("NCBIGene:1636", "biolink:same_as", "UniProtKB:P12821"),
                "attributes": {"infores:sri-node-normalizer": {"babel_relation": "gene_protein_conflation"}},
            },
            edge("RXCUI:861007", "biolink:has_active_ingredient", "CHEBI:6801"),
            edge("RXCUI:861007", "biolink:close_match", "CHEBI:6801"),
        ],
    )
    # Another source types CHEBI:6801 more coarsely and lists an id Babel doesn't know.
    refmet = _write_source(
        tmp_path,
        "refmet",
        [
            {
                "id": "RM:1",
                "categories": ["biolink:ChemicalEntity"],
                "provided_by": ["infores:refmet"],
                "equivalent_ids": ["RM:1", "CHEBI:6801"],
                "name": "Metformin",
            }
        ],
    )
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"babel": babel, "refmet": refmet},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    m = resolve_entities(config, biolink=None)

    assert m["NCBIGene:1636"] == m["HGNC:2707"] == m["UniProtKB:P12821"]
    assert m["RM:1"] == m["CHEBI:6801"]
    assert m["RXCUI:861007"] != m["CHEBI:6801"], "a drug/chemical relation merged a product into its ingredient"
    assert not (integrated / build_mod.BABEL_FACTS_FILENAME).exists()

    with jsonlines.open(config.integrated_nodes_path) as reader:
        by_id = {member: n for n in reader for member in n["equivalent_ids"]}
    assert "biolink:SmallMolecule" in by_id["CHEBI:6801"]["categories"]  # Babel's per-id category was used
    assert set(by_id["NCBIGene:1636"]["categories"]) == {"biolink:Gene", "biolink:Protein"}


def test_babel_evidence_regroups_stars_into_cliques_and_skips_drug_relations(tmp_path):
    import io

    from kraken.entity_resolution.build import _write_babel_evidence
    from kraken.entity_resolution.weights import ERWeights

    edges = [
        {"subject": "A:1", "predicate": "biolink:same_as", "object": "B:1"},
        {"subject": "A:1", "predicate": "biolink:same_as", "object": "C:1"},
        {
            "subject": "G:1",
            "predicate": "biolink:same_as",
            "object": "P:1",
            "attributes": {"infores:sri-node-normalizer": {"babel_relation": "gene_protein_conflation"}},
        },
        {"subject": "D:1", "predicate": "biolink:has_active_ingredient", "object": "A:1"},
        {"subject": "D:1", "predicate": "biolink:close_match", "object": "A:1"},
    ]
    path = tmp_path / "edges.jsonl"
    with jsonlines.open(path, "w") as writer:
        writer.write_all(edges)
    out = io.StringIO()
    _write_babel_evidence(path, ERWeights(), out)
    rows = [line.split("\t") for line in out.getvalue().splitlines()]
    pairs = {(a, b, kind) for a, b, _group, _weight, kind in rows}
    # the clique is whole again: B:1 and C:1 are linked to each other, not just to the hub
    assert pairs == {
        ("A:1", "B:1", "babel"),
        ("A:1", "C:1", "babel"),
        ("B:1", "C:1", "babel"),
        ("G:1", "P:1", "babel:gene_protein"),
    }
    assert {float(weight) for *_, weight, _kind in rows} == {ERWeights().equivalency_weight("babel")}
