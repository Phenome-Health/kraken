"""Hand-checked corrections of single wrong statements in a source's records (helpers/source_corrections.py).

The rule that matters: a correction touches only the named source's record. The same id from another source is a
separate claim and must come through whole.
"""

import json

import jsonlines
import pytest

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import NODE_DESCRIPTION, NODE_EQUIVALENT_IDS, NODE_SYNONYMS
from tests.helpers import stub_normalization


class _StubBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


class _Harmonizer(BaseHarmonizer):
    pass


def _harmonizer(source_infores: str) -> BaseHarmonizer:
    """Allocated without __init__, which reaches the network (Biolink toolkit, biomapper2)."""
    instance = object.__new__(_Harmonizer)
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


def test_the_same_record_from_another_source_is_left_alone():
    other = _harmonizer("infores:robokop-kg")
    assert "CHEBI:232965" in _econazole(other)[NODE_EQUIVALENT_IDS]
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
