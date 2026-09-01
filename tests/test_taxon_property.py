"""Tests for the single-valued `taxon` node property.

Biolink's taxon slot is single-valued, and KRAKEN's entity resolution enforces one taxon per merged entity, so
a node carries at most one. These go through `create_node` (the real path) rather than a helper, and stub the
Biolink client so the suite stays offline.
"""

import pytest

from kraken.harmonizers import base
from kraken.harmonizers.base import BaseHarmonizer
from kraken.schema import NodeModel
from kraken.utils.constants import NODE_TAXON


class _StubBiolink:
    """Only create_node's category filtering is needed; the real client downloads the Biolink model."""

    def filter_to_leaf_categories(self, categories):
        return list(categories)


class _Harmonizer(BaseHarmonizer):
    source_infores = "infores:test"


@pytest.fixture
def harmonizer() -> BaseHarmonizer:
    """Allocated without __init__, which builds a Biolink toolkit and a biomapper2 Normalizer (both reach the
    network). create_node needs only the stub client and the multi-taxon counters."""
    instance = object.__new__(_Harmonizer)
    instance.biolink = _StubBiolink()
    instance.multi_taxon_node_count = 0
    instance.multi_taxon_examples = []
    return instance


def _node(harmonizer, taxon, curie="NCBIGene:1"):
    return harmonizer.create_node(curie=curie, categories=["biolink:Gene"], provided_by="infores:test", taxon=taxon)


# ------------------------------- schema -------------------------------


def test_taxon_is_a_single_valued_string():
    """The property was list-valued (`taxa`) until entity resolution could guarantee one taxon per entity."""
    assert NODE_TAXON == "taxon"
    assert NodeModel.taxon.dtype is str
    assert NodeModel.taxon.inner_type is None
    assert not hasattr(NodeModel, "taxa")


# ------------------------------- one taxon -------------------------------


def test_single_taxon_is_stored_as_a_string(harmonizer):
    """However the source hands it over -- bare, listed, or as a set of repeats."""
    for supplied in ("NCBITaxon:9606", ["NCBITaxon:9606"], {"NCBITaxon:9606"}):
        assert _node(harmonizer, supplied)[NODE_TAXON] == "NCBITaxon:9606"
    assert harmonizer.multi_taxon_node_count == 0


def test_absent_taxon_omits_the_property(harmonizer):
    """Most nodes have no taxon at all -- valid, unlike having several."""
    for empty in (None, "", [], set(), [""], [None]):
        assert NODE_TAXON not in _node(harmonizer, empty)
    assert harmonizer.multi_taxon_node_count == 0


# ------------------------------- several taxa -------------------------------


def test_multiple_taxa_leave_the_node_untaxoned(harmonizer):
    """Several taxa means the source conflated distinct entities (translator labels the human gene CDH5 as
    human AND mouse). There's no way to tell which is meant, so the node gets none -- untaxoned nodes are
    wildcards for taxon merge guards, letting a source that knows the taxon supply it during resolution."""
    node = _node(harmonizer, ["NCBITaxon:9606", "NCBITaxon:10090"], curie="NCBIGene:1003")
    assert NODE_TAXON not in node
    assert harmonizer.multi_taxon_node_count == 1


def test_unsplit_delimited_taxon_value_leaves_the_node_untaxoned(harmonizer):
    """Translator supplies some taxa as one pipe-joined string alongside a second value. The node must end up
    untaxoned rather than storing the malformed CURIE as though it were a real taxon."""
    node = _node(harmonizer, ["NCBITaxon:9606|NCBITaxon:5476", "NCBITaxon:9606"], curie="NCBIGene:133")
    assert NODE_TAXON not in node
    assert harmonizer.multi_taxon_node_count == 1


def test_conflicts_are_counted_with_capped_examples(harmonizer):
    """Every occurrence is counted, but only a few are named, so a badly broken source can't flood the log."""
    for i in range(base.MAX_MULTI_TAXON_EXAMPLES + 3):
        _node(harmonizer, ["NCBITaxon:9606", "NCBITaxon:10090"], curie=f"NCBIGene:{i}")
    assert harmonizer.multi_taxon_node_count == base.MAX_MULTI_TAXON_EXAMPLES + 3
    assert len(harmonizer.multi_taxon_examples) == base.MAX_MULTI_TAXON_EXAMPLES


def test_example_names_the_node_and_its_taxa(harmonizer):
    """An example has to be enough to find the offending record back in the source."""
    _node(harmonizer, ["NCBITaxon:10116", "NCBITaxon:10090"], curie="NCBIGene:117063")
    example = harmonizer.multi_taxon_examples[0]
    assert "NCBIGene:117063" in example
    assert "NCBITaxon:10090" in example and "NCBITaxon:10116" in example
