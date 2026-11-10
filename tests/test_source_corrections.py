"""Hand-checked corrections of single wrong statements in a source's records (helpers/source_corrections.py).

The rule that matters: a correction touches only the named source's record. The same id from another source is a
separate claim and must come through whole.
"""

import json

import jsonlines
import pytest

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.kg2 import KG2Harmonizer
from kraken.utils.constants import NODE_DESCRIPTION, NODE_EQUIVALENT_IDS, NODE_SYNONYMS
from tests.helpers import stub_normalization


class _StubBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


class _Harmonizer(BaseHarmonizer):
    pass


def _harmonizer(source_infores: str, harmonizer_class: type[BaseHarmonizer] = _Harmonizer) -> BaseHarmonizer:
    """Allocated without __init__, which reaches the network (Biolink toolkit, biomapper2)."""
    instance = object.__new__(harmonizer_class)
    instance.source_infores = source_infores
    instance.biolink = _StubBiolink()
    instance.multi_taxon_node_count = 0
    instance.multi_taxon_examples = []
    instance.applied_corrections = set()
    instance.source_exclusions = set()
    instance.stripped_publications_count = 0
    instance.core_edge_props = {
        instance.subject_prop,
        instance.object_prop,
        instance.predicate_prop,
        instance.primary_ks_prop,
        instance.knowledge_level_prop,
        instance.agent_type_prop,
    }
    return stub_normalization(instance)


def _econazole(harmonizer: BaseHarmonizer) -> dict:
    """RefMet's econazole nitrate row, which lists N-succinyl-phenylalanine's ChEBI id (#29)."""
    return harmonizer.create_node(
        curie="RM:0048715",
        categories=["biolink:SmallMolecule"],
        provided_by=harmonizer.source_infores,
        name="Econazole nitrate",
        equivalent_ids=["RM:0048715", "CHEBI:232965", "PUBCHEM.COMPOUND:68589"],
    )


def _n_acetylasparagine(harmonizer: BaseHarmonizer) -> dict:
    """KG2's record of CHEBI:139582, which carries another drug's name (#29)."""
    return harmonizer.create_node(
        curie="CHEBI:139582",
        categories=["biolink:SmallMolecule"],
        provided_by=harmonizer.source_infores,
        name="N-alpha-acetyl-L-asparagine",
        synonyms=["N-acetylasparagine", "ACEGLUTAMIDE ALUMINIUM"],
    )


def test_a_correction_drops_the_wrong_statement_from_that_sources_record():
    refmet, kg2 = _harmonizer("infores:refmet"), _harmonizer("infores:rtx-kg2")
    assert sorted(_econazole(refmet)[NODE_EQUIVALENT_IDS]) == ["PUBCHEM.COMPOUND:68589", "RM:0048715"]
    assert sorted(_n_acetylasparagine(kg2)[NODE_SYNONYMS]) == ["N-acetylasparagine", "N-alpha-acetyl-L-asparagine"]
    assert refmet.applied_corrections == {("infores:refmet", "RM:0048715", "equivalent_id", "CHEBI:232965")}


def test_a_statement_the_source_no_longer_makes_is_reported_stale(caplog):
    """Once the source fixes its record the entry does nothing, and the log says to retire it."""
    refmet = _harmonizer("infores:refmet")
    refmet.create_node(
        curie="RM:0048715",
        categories=["biolink:SmallMolecule"],
        provided_by="infores:refmet",
        equivalent_ids=["RM:0048715", "PUBCHEM.COMPOUND:68589"],
    )
    refmet.log_source_corrections()
    assert "CHEBI:232965" in caplog.text and "retire" in caplog.text


def _cotinine_records(harmonizer: BaseHarmonizer, macrocycle_ids: list[str], cotinine_ids: list[str] | None) -> list:
    """A source's record of a screening-library macrocycle listing PathWhiz's cotinine, and its own cotinine (#29)."""
    records = [("CHEBI:127762", macrocycle_ids)] + ([("CHEBI:68641", cotinine_ids)] if cotinine_ids is not None else [])
    return [
        harmonizer.create_node(curie=curie, categories=["biolink:SmallMolecule"], provided_by="x", equivalent_ids=ids)
        for curie, ids in records
    ]


def test_a_moved_id_leaves_the_wrong_record_for_the_right_one(caplog):
    kg2 = _harmonizer("infores:rtx-kg2")
    macrocycle, cotinine = _cotinine_records(kg2, ["PathWhiz.Compound:824"], ["MESH:D003367"])
    assert "PathWhiz.Compound:824" not in macrocycle[NODE_EQUIVALENT_IDS]
    assert "PathWhiz.Compound:824" in cotinine[NODE_EQUIVALENT_IDS]
    kg2.log_source_corrections()
    assert "PathWhiz.Compound:824" not in caplog.text


def test_a_move_whose_right_record_is_gone_says_the_id_lost_its_record(caplog):
    """The move took the id off the wrong record, and nothing else gives it a record: it and its edges are gone."""
    kg2 = _harmonizer("infores:rtx-kg2")
    _cotinine_records(kg2, ["PathWhiz.Compound:824"], cotinine_ids=None)
    kg2.log_source_corrections()
    lost = [record.message for record in caplog.records if "no record to move to" in record.message]
    assert len(lost) == 1 and "PathWhiz.Compound:824" in lost[0]


def test_a_move_the_source_made_itself_is_retired_not_lost(caplog):
    """The wrong record no longer lists the id, so nothing was taken off it: the entry is stale, not a lost id."""
    kg2 = _harmonizer("infores:rtx-kg2")
    _cotinine_records(kg2, ["PUBCHEM.COMPOUND:54649749"], cotinine_ids=None)
    kg2.log_source_corrections()
    assert "no record to move to" not in caplog.text
    assert "PathWhiz.Compound:824" in caplog.text and "retire" in caplog.text


def test_a_right_record_that_already_lists_the_id_still_needs_the_move(caplog):
    """The wrong record still lists it, so the entry is not stale: retiring it would merge the two records."""
    kg2 = _harmonizer("infores:rtx-kg2")
    _cotinine_records(kg2, ["PathWhiz.Compound:824"], ["PathWhiz.Compound:824"])
    kg2.log_source_corrections()
    assert "PathWhiz.Compound:824" not in caplog.text


def _kg2_edge(*statements: str) -> dict:
    """A KG2 edge as KG2 writes it: its stored endpoints, and the original statements behind it."""
    return {
        "subject": "PathWhiz.Reaction:2021",
        "predicate": "biolink:has_input",
        "object": "CHEBI:127762",
        "primary_knowledge_source": "infores:pathwhiz",
        "knowledge_level": "knowledge_assertion",
        "agent_type": "manual_agent",
        "kg2_ids": list(statements),
    }


def _statement(subject: str, relation: str, object_: str, upstream: str) -> str:
    return "---".join([subject, relation, "None", "None", "None", object_, upstream])


PATHWHIZ_ON_PATHWHIZ_ID = _statement(
    "PathWhiz.Reaction:2021", "PathWhiz:has_left_element", "PathWhiz.Compound:824", "PathWhiz:"
)
PATHWHIZ_ON_CHEBI_ID = _statement("PathWhiz.Reaction:2021", "PathWhiz:has_left_element", "CHEBI:127762", "PathWhiz:")
CHEBI_ON_CHEBI_ID = _statement("CHEBI:127762", "rdfs:subClassOf", "CHEBI:24995", "OBO:chebi.owl")


def test_kg2_drops_the_original_statements_that_put_an_edge_on_the_wrong_node(tmp_path):
    """Integration places a KG2 edge by its original statements, so PathWhiz's cotinine edges keep only the one that
    names PathWhiz's own id, an edge with nothing else is dropped, and other sources' statements on the id stay."""
    edges_in, edges_out = tmp_path / "edges.jsonl", tmp_path / "out.jsonl"
    edges = [
        _kg2_edge(PATHWHIZ_ON_PATHWHIZ_ID, PATHWHIZ_ON_CHEBI_ID),
        _kg2_edge(PATHWHIZ_ON_CHEBI_ID),
        {**_kg2_edge(CHEBI_ON_CHEBI_ID), "primary_knowledge_source": "infores:chebi"},
    ]
    edges_in.write_text("\n".join(json.dumps(edge) for edge in edges))
    kg2 = _harmonizer("infores:rtx-kg2", KG2Harmonizer)
    kg2._harmonize_edges(edges_in, edges_out)
    with jsonlines.open(edges_out) as reader:
        statements = [edge["attributes"]["infores:rtx-kg2"]["kg2pre_ids"] for edge in reader]
    assert statements == [[PATHWHIZ_ON_PATHWHIZ_ID], [CHEBI_ON_CHEBI_ID]]
    assert ("infores:rtx-kg2", "CHEBI:127762", "PathWhiz:") in kg2.applied_corrections  # else reported stale


def test_the_same_record_from_another_source_is_left_alone():
    other = _harmonizer("infores:robokop-kg")
    assert "CHEBI:232965" in _econazole(other)[NODE_EQUIVALENT_IDS]
    macrocycle, cotinine = _cotinine_records(other, ["PathWhiz.Compound:824"], ["MESH:D003367"])
    assert "PathWhiz.Compound:824" in macrocycle[NODE_EQUIVALENT_IDS]
    assert "PathWhiz.Compound:824" not in cotinine[NODE_EQUIVALENT_IDS]
    assert "ACEGLUTAMIDE ALUMINIUM" in _n_acetylasparagine(other)[NODE_SYNONYMS]


@pytest.mark.parametrize("source, wrong_edge_is_dropped", [("infores:rtx-kg2", True), ("infores:robokop-kg", False)])
def test_an_edge_correction_drops_only_that_sources_edge(tmp_path, source, wrong_edge_is_dropped):
    """Reactome's phenylacetate same_as the ester phenyl acetate, as KG2 republishes it (#29)."""
    edge = {
        "predicate": "biolink:same_as",
        "primary_knowledge_source": "infores:reactome",
        "knowledge_level": "knowledge_assertion",
        "agent_type": "manual_agent",
    }
    edges_in = tmp_path / "edges.jsonl"
    edges_in.write_text(
        "\n".join(
            json.dumps({**edge, "subject": "REACT:R-ALL-177153", "object": object_})
            for object_ in ("CHEBI:8082", "CHEBI:14779")
        )
    )
    edges_out = tmp_path / "out.jsonl"
    _harmonizer(source)._harmonize_edges(edges_in, edges_out)
    with jsonlines.open(edges_out) as reader:
        objects = {edge["object"] for edge in reader}
    assert objects == ({"CHEBI:14779"} if wrong_edge_is_dropped else {"CHEBI:8082", "CHEBI:14779"})


@pytest.mark.parametrize(
    "description, dropped",
    [
        (
            "(alpha-D-mannosyl)7-beta-D-mannosyl-diacetylchitobiosyl-L-asparagine, isoform B (protein) is an anisole.",
            True,
        ),
        ("4-Iodoanisole is an iodo-substituted anisole.", False),
    ],
    ids=["the audited text", "a repaired text"],
)
def test_a_description_goes_only_while_it_says_what_was_audited(description, dropped):
    kg2 = _harmonizer("infores:rtx-kg2")
    node = kg2.create_node(
        curie="CHEBI:174192", categories=["biolink:SmallMolecule"], provided_by="x", description=description
    )
    assert (NODE_DESCRIPTION not in node) == dropped
