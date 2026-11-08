"""An id UMLS minted itself joins its concept's node, without disturbing anything else (see umls_minted_codes)."""

from pathlib import Path
from types import SimpleNamespace

import jsonlines
import pytest

from kraken.entity_resolution.build import resolve_entities
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo
from kraken.entity_resolution.materialize import PrefixRanking
from kraken.entity_resolution.umls_minted_codes import attach_orphan_codes

CLIQUE = ["MONDO:0005148", "UMLS:C0011860", "DOID:9352"]


def _kg2_node(curie: str, category: str, equivalent_ids: list[str] | None = None) -> dict:
    return {
        "id": curie,
        "name": curie,
        "categories": [category],
        "equivalent_ids": equivalent_ids or [curie],
        "provided_by": ["infores:rtx-kg2"],
    }


def _xref_edge(code: str, concept: str) -> dict:
    """kg2 stores it on the canonical id, and keeps the original endpoints in kg2pre_ids."""
    return {
        "subject": code,
        "object": CLIQUE[0],
        "predicate": "biolink:close_match",
        "primary_knowledge_source": "infores:omim",
        "attributes": {
            "infores:rtx-kg2": {
                "kg2pre_ids": [f"{code}---UMLS:xref---None---None---None---{concept}---umls_source:OMIM"]
            }
        },
    }


def _resolve(
    tmp_path: Path,
    kg2_nodes: list[dict],
    kg2_edges: list[dict],
    clique: list[str] = CLIQUE,
    categories: tuple[str, ...] = ("biolink:Disease",),
) -> dict[str, str]:
    harmonized = {}
    babel_nodes = [{"id": c, "name": c, "categories": list(categories), "provided_by": ["x"]} for c in clique]
    babel_edges = [
        {"subject": clique[0], "predicate": "biolink:same_as", "object": c, "primary_knowledge_source": "x"}
        for c in clique[1:]
    ]
    for source, nodes, edges in (("babel", babel_nodes, babel_edges), ("kg2", kg2_nodes, kg2_edges)):
        directory = tmp_path / "harmonized" / source
        directory.mkdir(parents=True)
        with jsonlines.open(directory / "nodes.jsonl", "w") as writer:
            writer.write_all(nodes)
        with jsonlines.open(directory / "edges.jsonl", "w") as writer:
            writer.write_all(edges)
        harmonized[source] = (directory / "nodes.jsonl", directory / "edges.jsonl")
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    return resolve_entities(config, biolink=None)


def test_code_with_no_other_tie_joins_its_concept(tmp_path):
    """The type 2 diabetes split: its only evidence is one close_match, below tau."""
    code = "OMIM:MTHU069260"
    mapping = _resolve(tmp_path, [_kg2_node(code, "biolink:Disease")], [_xref_edge(code, "UMLS:C0011860")])
    assert mapping[code] == mapping["MONDO:0005148"]


def test_code_joins_when_no_pair_reached_tau(tmp_path):
    """No pair means nothing to cluster, not nothing to attach. The default ranking, which a Phenomenon gets, puts
    UMLS above OMIM, so the code keeps the concept's id."""
    code, concept = "OMIM:MTHU000005", "UMLS:C0000005"
    nodes, edges = [_kg2_node(code, "biolink:Phenomenon")], [_xref_edge(code, concept)]
    mapping = _resolve(tmp_path, nodes, edges, clique=[concept], categories=("biolink:Phenomenon",))
    assert mapping[code] == mapping[concept]


def test_guardrails_can_refuse_it(tmp_path):
    code = "OMIM:MTHU000001"
    mapping = _resolve(tmp_path, [_kg2_node(code, "biolink:Gene")], [_xref_edge(code, "UMLS:C0011860")])
    assert mapping[code] != mapping["MONDO:0005148"]


@pytest.mark.parametrize(
    ("clique", "category", "taxon"),
    [
        pytest.param(["HGNC:1100", "UMLS:C0376571"], "biolink:Gene", "NCBITaxon:9606", id="the node's own taxon"),
        pytest.param(["HGNC:1100", "UMLS:C0376571"], "biolink:Gene", "NCBITaxon:10090", id="another species"),
        pytest.param(CLIQUE, "biolink:Disease", "NCBITaxon:9606", id="a node with no taxon"),
    ],
)
def test_code_with_a_taxon_stays_apart(tmp_path, clique, category, taxon):
    code = "OMIM:MTHU000006"
    nodes, edges = [{**_kg2_node(code, category), "taxon": taxon}], [_xref_edge(code, clique[1])]
    mapping = _resolve(tmp_path, nodes, edges, clique=clique, categories=(category,))
    assert mapping[code] != mapping[clique[0]]


def test_the_guardrails_decide_too():
    """Nothing else here would refuse a second id of an enforced prefix."""
    code, concept = "OMIM:MTHU000009", "UMLS:C0000009"
    curie_to_cluster = {concept: 0, "OMIM:600001": 0}
    attached, _ = attach_orphan_codes(
        {code: concept},
        curie_to_cluster,
        info_of=lambda curie: NodeInfo(curie=curie, branches=frozenset({"disease_pheno"})),
        categories_of={curie: ("biolink:Disease",) for curie in (code, concept, "OMIM:600001")},
        label_of=lambda curie: curie,
        ranking=PrefixRanking.load(),
        guardrail_config=GuardrailConfig(enforced_prefixes=frozenset({"OMIM"})),
        next_cluster_id=1,
    )
    assert attached == set()


def test_code_does_not_join_a_chemical(tmp_path):
    """Everything else would let it in: CHEBI keeps the lead and the types match. But a chemical's name can be a
    spelling its members vote on, and the code's name would be one more vote."""
    code, clique = "OMIM:MTHU000007", ["CHEBI:15365", "UMLS:C0004057"]
    nodes, edges = [_kg2_node(code, "biolink:SmallMolecule")], [_xref_edge(code, clique[1])]
    mapping = _resolve(tmp_path, nodes, edges, clique=clique, categories=("biolink:SmallMolecule",))
    assert mapping[code] != mapping["CHEBI:15365"]


def test_code_that_would_change_the_node_stays_apart(tmp_path):
    """It may add only itself: OMIM ranks above HP, and kg2 types a code from its concept's UMLS semantic type."""
    code = "OMIM:MTHU000093"  # 'Abnormal pinnae'
    nodes, edges = [_kg2_node(code, "biolink:PhenotypicFeature")], [_xref_edge(code, "UMLS:C1846460")]
    hp_node = ["HP:0000356", "UMLS:C1846460"]
    would_lead = _resolve(tmp_path / "id", nodes, edges, clique=hp_node, categories=("biolink:PhenotypicFeature",))
    assert would_lead[code] != would_lead["HP:0000356"]

    code = "OMIM:MTHU000003"
    nodes, edges = [_kg2_node(code, "biolink:PhenotypicFeature")], [_xref_edge(code, "UMLS:C0011860")]
    would_retype = _resolve(tmp_path / "type", nodes, edges)
    assert would_retype[code] != would_retype["MONDO:0005148"]


def test_code_that_would_hand_the_node_to_another_id_stays_apart(tmp_path):
    """A ClinicalFinding also sits in the clinical family, which ranks UMLS above DOID. A Disease-only code would narrow
    the node to the disease family, where DOID ranks first, and so change its id and name without leading it."""
    code, clique = "OMIM:MTHU000004", ["DOID:9352", "UMLS:C0011860"]
    nodes, edges = [_kg2_node(code, "biolink:Disease")], [_xref_edge(code, "UMLS:C0011860")]
    mapping = _resolve(tmp_path, nodes, edges, clique=clique, categories=("biolink:Disease", "biolink:ClinicalFinding"))
    assert mapping["DOID:9352"] == "UMLS:C0011860"
    assert mapping[code] != mapping["DOID:9352"]


def test_code_that_would_rename_the_node_stays_apart(tmp_path):
    """A node whose top id has no Babel label is named after its best-ranked member with a name, which the code would
    become: OMIM:MTHU069324 would have renamed the 2.4.0 node of OMIM:MTHU069320, now named after UMLS:C3281069."""
    leader, concept, code = "OMIM:MTHU069320", "UMLS:C3281069", "OMIM:MTHU069324"
    nodes = [_kg2_node(concept, "biolink:Disease", [concept, leader]), _kg2_node(code, "biolink:Disease")]
    _resolve(tmp_path, nodes, [_xref_edge(code, concept)], clique=[concept])
    node = next(n for n in jsonlines.open(tmp_path / "integrated" / "nodes.jsonl") if n["id"] == leader)
    assert node["name"] == concept  # every fixture node is named after its id


def test_the_node_keeps_its_description(tmp_path):
    """MONDO and DOID have none here, so the node takes the best-ranked member's that has one: UMLS's, and not the
    joining code's semantic type, though OMIM ranks above UMLS."""
    code, concept = "OMIM:MTHU069260", "UMLS:C0011860"
    nodes = [
        {**_kg2_node(concept, "biolink:Disease"), "description": "A type of diabetes mellitus."},
        {**_kg2_node(code, "biolink:Disease"), "description": "UMLS Semantic Type: STY:T047"},
    ]
    mapping = _resolve(tmp_path, nodes, [_xref_edge(code, concept)])
    assert mapping[code] == mapping["MONDO:0005148"]
    node = next(n for n in jsonlines.open(tmp_path / "integrated" / "nodes.jsonl") if n["id"] == "MONDO:0005148")
    assert node["description"] == "A type of diabetes mellitus."


def test_a_cluster_opened_for_the_concept_is_not_handed_out_again():
    """The concept was in no cluster, so the code opened one. The next cluster made must not get its id."""
    code, concept = "OMIM:MTHU000026", "MONDO:0000002"
    curie_to_cluster = {"MONDO:0000001": 0}
    attached, next_cluster_id = attach_orphan_codes(
        {code: concept},
        curie_to_cluster,
        info_of=lambda curie: NodeInfo(curie=curie, branches=frozenset({"disease_pheno"})),
        categories_of={code: ("biolink:Disease",), concept: ("biolink:Disease",)},
        label_of=lambda curie: None if curie == code else curie,
        ranking=PrefixRanking.load(),
        guardrail_config=GuardrailConfig(),
        next_cluster_id=1,
    )
    assert attached == {code}
    assert curie_to_cluster[code] == curie_to_cluster[concept] == 1
    assert next_cluster_id == 2


def test_code_that_already_shares_a_node_stays_there(tmp_path):
    """Its other tie disagrees with the concept. Attaching it would pull one cluster's id into another."""
    code, sibling = "OMIM:MTHU000002", "CHV:0000000002"
    nodes = [_kg2_node(sibling, "biolink:Disease", [sibling, code])]
    mapping = _resolve(tmp_path, nodes, [_xref_edge(code, "UMLS:C0011860")])
    assert mapping[code] == mapping[sibling]
    assert mapping[code] != mapping["MONDO:0005148"]


def test_other_vocabularies_xrefs_stay_weak(tmp_path):
    """Other vocabularies' codes can sit in several concepts: a code in two would fuse them."""
    code = "CHV:0000000546"
    mapping = _resolve(tmp_path, [_kg2_node(code, "biolink:Disease")], [_xref_edge(code, "UMLS:C0011860")])
    assert mapping[code] != mapping["MONDO:0005148"]
