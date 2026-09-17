"""Tests for the Babel harmonizer, run over a miniature release directory.

The harmonizer is allocated without __init__ (which builds a Biolink toolkit and a biomapper2 Normalizer, both
of which reach the network); see ``_harmonizer`` for the state it needs instead.
"""

import json
from collections import Counter, defaultdict

import jsonlines

from kraken.harmonizers.babel import BabelHarmonizer
from tests.helpers import stub_normalization, write_test_taxdump


class _LeafBiolink:
    """Every category is already a leaf -- enough for create_node."""

    def filter_to_leaf_categories(self, categories):
        return list(categories) if isinstance(categories, (list, set, tuple)) else [categories]


def _identifier(curie, label="", taxa=(), descriptions=()):
    return {"i": curie, "l": label, "d": list(descriptions), "t": list(taxa)}


def _clique(category, identifiers, taxa=()):
    # Key order matches Babel's files: the harmonizer reads "taxa" off the end of the line.
    return {"type": category, "ic": None, "identifiers": identifiers, "preferred_name": "", "taxa": list(taxa)}


def _write_release(tmp_path, compendia, gene_protein=(), drug_chemical=(), relations=()):
    release = tmp_path / "release"
    (release / "compendia").mkdir(parents=True)
    (release / "conflation").mkdir()
    for name, cliques in compendia.items():
        (release / "compendia" / f"{name}.txt").write_text("".join(json.dumps(c) + "\n" for c in cliques))
    (release / "conflation" / "GeneProtein.txt").write_text("".join(json.dumps(g) + "\n" for g in gene_protein))
    (release / "conflation" / "DrugChemical.txt").write_text("".join(json.dumps(g) + "\n" for g in drug_chemical))
    (release / "conflation" / "DrugChemical_concords.tsv").write_text(
        "subj\tpred\tobj\n" + "".join("\t".join(row) + "\n" for row in relations)
    )
    # 9606 is a species (in the allowlist); 12345 is a human strain-level child that must roll up to it;
    # 999 is a species that is not allowlisted.
    write_test_taxdump(
        release,
        nodes=[("1", "1", "no rank"), ("9606", "1", "species"), ("12345", "9606", "strain"), ("999", "1", "species")],
        names={"9606": "Homo sapiens", "12345": "human strain", "999": "Nobody cares"},
    )
    return release


def _harmonizer() -> BabelHarmonizer:
    harmonizer = object.__new__(BabelHarmonizer)
    harmonizer.source_infores = "infores:sri-node-normalizer"
    harmonizer.biolink = _LeafBiolink()
    harmonizer.name_override_count = 0
    harmonizer.multi_taxon_node_count = 0
    harmonizer.multi_taxon_examples = []
    stub_normalization(harmonizer)
    harmonizer.taxonomy = None
    harmonizer._taxon_allowed = {}
    harmonizer._species_curie = {}
    harmonizer.gene_leaders = set()
    harmonizer.protein_leaders = set()
    harmonizer.drug_chemical_ids = set()
    harmonizer.drug_chemical_leader = {}
    harmonizer.stats = defaultdict(Counter)
    return harmonizer


def _run(tmp_path, **release):
    release_dir = _write_release(tmp_path, **release)
    harmonizer = _harmonizer()
    nodes_path, edges_path = tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl"
    harmonizer.harmonize(nodes_path, edges_path, input_file=release_dir)
    with jsonlines.open(nodes_path) as nodes, jsonlines.open(edges_path) as edges:
        return harmonizer, {n["id"]: n for n in nodes}, list(edges)


def _triples(edges):
    return {(e["subject"], e["predicate"], e["object"]) for e in edges}


def test_every_identifier_is_a_node_with_its_own_facts_and_cliques_are_same_as_stars(tmp_path):
    _, nodes, edges = _run(
        tmp_path,
        compendia={
            "Disease": [
                _clique(
                    "biolink:Disease",
                    [
                        _identifier("MONDO:1", "renal tubular dysgenesis", descriptions=["A kidney disorder."]),
                        _identifier("DOID:2", "RTD"),
                        _identifier("UMLS:C3"),
                    ],
                )
            ]
        },
    )
    assert set(nodes) == {"MONDO:1", "DOID:2", "UMLS:C3"}
    assert nodes["DOID:2"]["name"] == "RTD"  # the id's own label, not the clique's preferred name
    assert nodes["MONDO:1"]["description"] == "A kidney disorder."
    assert "name" not in nodes["UMLS:C3"]  # no label of its own -> none borrowed
    assert all(n["categories"] == ["biolink:Disease"] for n in nodes.values())
    assert all(n["equivalent_ids"] == [n["id"]] for n in nodes.values())  # the clique lives in the edges
    assert _triples(edges) == {("MONDO:1", "biolink:same_as", "DOID:2"), ("MONDO:1", "biolink:same_as", "UMLS:C3")}


def test_gene_protein_cliques_are_scoped_by_taxon_and_conflated_by_same_as(tmp_path):
    harmonizer, nodes, edges = _run(
        tmp_path,
        compendia={
            "Gene": [
                _clique(
                    "biolink:Gene",
                    [_identifier("NCBIGene:1636", "ACE", ["NCBITaxon:9606"]), _identifier("HGNC:2707", "ACE")],
                    taxa=["NCBITaxon:9606"],
                ),
                _clique("biolink:Gene", [_identifier("NCBIGene:5", "x", ["NCBITaxon:999"])], taxa=["NCBITaxon:999"]),
                _clique("biolink:Gene", [_identifier("UMLS:C9", "a gene concept")]),  # no taxon: kept
            ],
            "Protein": [
                _clique(
                    "biolink:Protein",
                    [_identifier("UniProtKB:P12821", "ACE", ["NCBITaxon:12345"])],
                    taxa=["NCBITaxon:12345"],
                ),
            ],
        },
        gene_protein=[["NCBIGene:1636", "UniProtKB:P12821"], ["NCBIGene:5", "UniProtKB:Q0"]],
    )
    assert "NCBIGene:5" not in nodes  # its organism isn't allowlisted
    assert "UMLS:C9" in nodes
    assert nodes["UniProtKB:P12821"]["taxon"] == "NCBITaxon:9606"  # strain rolled up to species
    conflation = [
        e
        for e in edges
        if e.get("attributes", {}).get("infores:sri-node-normalizer", {}).get("babel_relation")
        == "gene_protein_conflation"
    ]
    assert _triples(conflation) == {("NCBIGene:1636", "biolink:same_as", "UniProtKB:P12821")}  # not NCBIGene:5's
    assert harmonizer.stats["Gene"]["cliques_dropped_taxon"] == 1


def test_structure_only_chemical_cliques_are_dropped(tmp_path):
    _, nodes, _ = _run(
        tmp_path,
        compendia={
            "SmallMolecule": [
                _clique(
                    "biolink:SmallMolecule",
                    [_identifier("PUBCHEM.COMPOUND:1", "some structure"), _identifier("INCHIKEY:AAA")],
                ),
                # CAS and ChEMBL count as curation: these cliques are the combination products and assay
                # compounds the aggregators conflate, so Babel must know them to keep them apart
                _clique(
                    "biolink:SmallMolecule",
                    [_identifier("PUBCHEM.COMPOUND:2", "a registered structure"), _identifier("CAS:338392-03-3")],
                ),
                _clique(
                    "biolink:SmallMolecule",
                    [
                        _identifier("PUBCHEM.COMPOUND:46861711", "Linagliptin; METformin Hydrochloride"),
                        _identifier("CHEMBL.COMPOUND:CHEMBL4879149"),
                    ],
                ),
                _clique(
                    "biolink:SmallMolecule",
                    [_identifier("CHEBI:6801", "metformin"), _identifier("PUBCHEM.COMPOUND:4091")],
                ),
            ]
        },
    )
    assert set(nodes) == {
        "CHEBI:6801",
        "PUBCHEM.COMPOUND:4091",
        "PUBCHEM.COMPOUND:2",
        "CAS:338392-03-3",
        "PUBCHEM.COMPOUND:46861711",
        "CHEMBL.COMPOUND:CHEMBL4879149",
    }


def test_drug_chemical_relations_are_typed_edges_with_close_match_only_where_no_relation_connects(tmp_path):
    _, _, edges = _run(
        tmp_path,
        compendia={
            "SmallMolecule": [
                _clique("biolink:SmallMolecule", [_identifier("CHEBI:6801", "metformin"), _identifier("RXCUI:6809")]),
            ],
            "ChemicalEntity": [_clique("biolink:ChemicalEntity", [_identifier("RXCUI:235743", "metformin HCl")])],
            "Drug": [
                _clique("biolink:Drug", [_identifier("RXCUI:861007", "metformin HCl 500 MG Oral Tablet")]),
                _clique("biolink:Drug", [_identifier("RXCUI:7", "Glucophage")]),
            ],
        },
        # The group lists clique preferred ids; the relations use whichever id RxNorm has (RXCUI:6809 is not
        # CHEBI:6801's preferred id, but is in its clique).
        drug_chemical=[["CHEBI:6801", "RXCUI:235743", "RXCUI:861007", "RXCUI:7"]],
        relations=[
            ("RXCUI:861007", "has_precise_active_ingredient", "RXCUI:235743"),
            ("RXCUI:235743", "has_form", "RXCUI:6809"),
            ("RXCUI:861007", "has_ingredient", "RXCUI:404"),  # endpoint in no kept clique -> skipped
        ],
    )
    typed = [e for e in edges if e["predicate"] != "biolink:same_as"]
    assert _triples(typed) == {
        ("RXCUI:861007", "biolink:has_active_ingredient", "RXCUI:235743"),
        ("RXCUI:235743", "biolink:close_match", "RXCUI:6809"),
        # the only member no relation reaches from the group's preferred id
        ("CHEBI:6801", "biolink:close_match", "RXCUI:7"),
    }
    relation_of = {(e["subject"], e["object"]): e["attributes"]["infores:sri-node-normalizer"] for e in typed}
    assert relation_of[("RXCUI:861007", "RXCUI:235743")] == {"babel_relation": "has_precise_active_ingredient"}
    assert relation_of[("CHEBI:6801", "RXCUI:7")] == {"babel_relation": "drug_chemical_conflation"}


def test_unchanged_ids_are_not_cached_but_failures_are_reported_once(tmp_path):
    harmonizer = _harmonizer()

    class _RejectsFoo:
        def get_curies(self, local_ids_dict, **_kwargs):
            ((vocab, local_id),) = local_ids_dict.items()
            return ({}, {}, {vocab}) if vocab == "FOO" else ({f"{vocab}:{local_id}": ""}, {}, set())

    harmonizer.normalizer = _RejectsFoo()
    for _ in range(3):
        assert harmonizer.normalize_curie("MONDO:1") == "MONDO:1"
        assert harmonizer.normalize_curie("FOO:1") == "FOO:1"
    assert "MONDO:1" not in harmonizer.normalized_id_map
    assert harmonizer.unrecognized_vocab_prefixes["FOO"]["count"] == 1


def test_a_lone_identifier_with_no_label_and_no_taxon_is_dropped(tmp_path):
    _, nodes, _ = _run(
        tmp_path,
        compendia={
            "Gene": [
                _clique("biolink:Gene", [_identifier("ENSEMBL:ENSCJAG00000078489")]),  # nothing but an id
                _clique("biolink:Gene", [_identifier("MGI:8135641", "Rr673509")]),  # labeled: kept
            ],
            "Drug": [
                _clique("biolink:Drug", [_identifier("UMLS:C1243498")]),
                # unlabeled members are fine inside a real clique
                _clique("biolink:Drug", [_identifier("RXCUI:1726214"), _identifier("UMLS:C1618329")]),
            ],
        },
    )
    assert set(nodes) == {"MGI:8135641", "RXCUI:1726214", "UMLS:C1618329"}
