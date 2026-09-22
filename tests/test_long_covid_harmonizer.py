"""long-covid: every source column survives, and a symptom several tables report is one edge carrying all of them."""

import json

from kraken.harmonizers.long_covid import (
    COMPOUND_FILE,
    PROTEIN_FILE,
    LongCovidHarmonizer,
)
from tests.helpers import stub_normalization

SOURCE_ID = "isb-long-covid"

INCOV_SYMPTOMS = """disease_id,symptom_name,symptom_id,percentage_total_diagnosed_at_time_point1
DOID:0080848,Fatigue,D005221,0.88
,,,
,,,`
"""
INSIGHT_SYMPTOMS = """disease_id,symptom_name,symptom_id,incidence_rate_of_new_symptoms_among_OneFlorida_patients
DOID:0080848,Fatigue,D005221,0.05
"""
RECOVER_SYMPTOMS = """disease_id,symptom_name,symptom_id,symptom_frequency_among_overall_recover_patients
DOID:0080848,Cough,D003371,0.33
DOID:0080848,Urinary Bladder,D001743,0.3
"""
PROTEINS = """disease_id,gene_name,protein_name,protein_id,edge_type,pmid_list,preprint_list
DOID:0080848,CCL19,C-C motif chemokine 19,Q99731,INCREASEDIN_PiD,37748514,
DOID:0080848,CCL19,C-C motif chemokine 19,Q99731,INCREASEDIN_PiD,36947108,
DOID:0080848,IL5,Interleukin-5,P05113,DECREASEDIN_PdD,37748514,
"""
COMPOUNDS = """disease_id,compound_id,percentage_total_diagnosed_at_time_point1
DOID:0080848,CHEMBL196,0.04
,,`
"""


class _LeafBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


def _harmonizer() -> LongCovidHarmonizer:
    """A harmonizer allocated without __init__ (which builds a Biolink toolkit and biomapper2 Normalizer)."""
    harmonizer = object.__new__(LongCovidHarmonizer)
    harmonizer.source_infores = SOURCE_ID
    harmonizer.biolink = _LeafBiolink()
    harmonizer.name_override_count = 0
    harmonizer.multi_taxon_node_count = 0
    harmonizer.multi_taxon_examples = []
    return stub_normalization(harmonizer)


def _run(tmp_path) -> tuple[dict, dict]:
    for filename, content in {
        "INCOV_long_covid_symptom_edge.csv": INCOV_SYMPTOMS,
        "INSIGHT_OneFlorida_long_covid_symptom_edge.csv": INSIGHT_SYMPTOMS,
        "RECOVER_long_covid_symptom_edge.csv": RECOVER_SYMPTOMS,
        PROTEIN_FILE: PROTEINS,
        COMPOUND_FILE: COMPOUNDS,
    }.items():
        (tmp_path / filename).write_text(content)
    nodes_path, edges_path = tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl"
    _harmonizer().harmonize(nodes_path, edges_path, input_file=tmp_path)
    nodes = {n["id"]: n for n in map(json.loads, nodes_path.read_text().splitlines())}
    edge_list = map(json.loads, edges_path.read_text().splitlines())
    edges = {(e["subject"], e["predicate"], e["object"]): e for e in edge_list}
    return nodes, edges


def test_symptom_reported_by_two_tables_is_one_edge_with_both_tables_data(tmp_path):
    _, edges = _run(tmp_path)
    edge = edges[("DOID:0080848", "biolink:has_phenotype", "MESH:D005221")]
    attributes = edge["attributes"][SOURCE_ID]
    assert attributes["percentage_total_diagnosed_at_time_point1"] == 0.88  # INCOV
    assert attributes["incidence_rate_of_new_symptoms_among_OneFlorida_patients"] == 0.05  # INSIGHT/OneFlorida
    assert len(attributes["source_files"]) == 2
    assert edge["publications"] == ["PMID:35216672", "PMID:37029117", "PMID:36785842"]
    assert (edge["knowledge_level"], edge["agent_type"]) == ("statistical_association", "data_analysis_pipeline")


def test_protein_edge_type_selects_correlation_direction_and_pmids_merge(tmp_path):
    nodes, edges = _run(tmp_path)
    ccl19 = edges[("DOID:0080848", "biolink:positively_correlated_with", "UniProtKB:Q99731")]
    assert ccl19["publications"] == ["PMID:37748514", "PMID:36947108"]
    assert ccl19["attributes"][SOURCE_ID]["edge_type"] == "INCREASEDIN_PiD"
    assert ("DOID:0080848", "biolink:negatively_correlated_with", "UniProtKB:P05113") in edges
    assert nodes["UniProtKB:Q99731"]["name"] == "C-C motif chemokine 19"
    assert "CCL19" in nodes["UniProtKB:Q99731"]["synonyms"]


def test_compound_edge_points_at_disease_and_blank_rows_are_skipped(tmp_path):
    nodes, edges = _run(tmp_path)
    edge = edges[("CHEMBL.COMPOUND:CHEMBL196", "biolink:applied_to_treat", "DOID:0080848")]
    assert edge["attributes"][SOURCE_ID]["percentage_total_diagnosed_at_time_point1"] == 0.04
    assert len(edges) == 6  # 3 symptoms + 2 proteins + 1 compound; the trailing blank rows add nothing
    assert len(nodes) == 7


def test_non_phenotype_mesh_terms_are_named_things(tmp_path):
    nodes, _ = _run(tmp_path)
    assert nodes["MESH:D003371"]["categories"] == ["biolink:DiseaseOrPhenotypicFeature"]  # Cough
    assert nodes["MESH:D001743"]["categories"] == ["biolink:NamedThing"]  # Urinary Bladder
