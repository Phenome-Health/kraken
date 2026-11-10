"""KG2's UMLS semantic types leave the description and become an attribute (#17). Examples are real KG2 shapes."""

import pytest

from kraken.harmonizers.kg2 import KG2Harmonizer, split_umls_semantic_types
from kraken.utils.constants import NODE_ATTRIBUTES, NODE_DESCRIPTION
from tests.helpers import stub_normalization


@pytest.mark.parametrize(
    "description, text, semantic_types",
    [
        ("UMLS Semantic Type: STY:T121", "", ["STY:T121"]),
        ("UMLS Semantic Type: STY:T121; UMLS Semantic Type: STY:T109", "", ["STY:T109", "STY:T121"]),
        ("UMLS Semantic Type: STY:T204; ; ;", "", ["STY:T204"]),
        ("An agent. (NCI04); UMLS Semantic Type: STY:T200", "An agent. (NCI04)", ["STY:T200"]),
        (
            "UMLS Semantic Type: STY:T109; UMLS Semantic Type: STY:T121; Domiphen bromide is an antiseptic.",
            "Domiphen bromide is an antiseptic.",
            ["STY:T109", "STY:T121"],
        ),
        (
            "In the bone marrow.; UMLS Semantic Type: STY:T033; UMLS Semantic Type: STY:T047; Increased count.;",
            "In the bone marrow.; Increased count.",
            ["STY:T033", "STY:T047"],
        ),
        # an HTML entity's ";" is not a separator (HP:0006379)
        (
            "UMLS Semantic Type: STY:T033; New term from UBERON: &apos;proximal tibia&apos;",
            "New term from UBERON: &apos;proximal tibia&apos;",
            ["STY:T033"],
        ),
        ("A &lt;b&gt; tag.&#39;; UMLS Semantic Type: STY:T047", "A &lt;b&gt; tag.&#39;", ["STY:T047"]),
    ],
)
def test_semantic_types_leave_the_description(description, text, semantic_types):
    assert split_umls_semantic_types(description) == (text, semantic_types)


def test_the_sources_own_empty_parts_stay():
    """Only empty parts touching a removed type or an end go; the rest is the source's punctuation."""
    description = "A sudden onset of abdominal pain.; ; Examples: one abdomen.; UMLS Semantic Type: STY:T029"
    assert split_umls_semantic_types(description) == (
        "A sudden onset of abdominal pain.; ; Examples: one abdomen.",
        ["STY:T029"],
    )
    assert split_umls_semantic_types("No types here.; ; Still none.") == ("No types here.; ; Still none.", [])


class _StubBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


@pytest.fixture
def kg2() -> KG2Harmonizer:
    """Allocated without __init__, which reaches the network (Biolink toolkit, biomapper2)."""
    instance = object.__new__(KG2Harmonizer)
    instance.source_infores = "infores:rtx-kg2"
    instance.biolink = _StubBiolink()
    instance.multi_taxon_node_count = 0
    instance.multi_taxon_examples = []
    instance.applied_corrections = set()
    instance.name_override_count = 0
    instance.core_node_props = {"id", "all_categories", "equivalent_curies", "name", "description", "iri", "all_names"}
    return stub_normalization(instance)


def _kg2_node(kg2: KG2Harmonizer, description: str) -> dict:
    return kg2._harmonize_node({"id": "UMLS:C1", "all_categories": ["biolink:Disease"], "description": description})


def test_kg2_moves_the_semantic_types_to_an_attribute(kg2):
    node = _kg2_node(kg2, "UMLS Semantic Type: STY:T121; UMLS Semantic Type: STY:T109")
    assert NODE_DESCRIPTION not in node
    assert node[NODE_ATTRIBUTES]["infores:rtx-kg2"]["umls_semantic_types"] == ["STY:T109", "STY:T121"]


def test_kg2_keeps_the_real_text(kg2):
    node = _kg2_node(kg2, "An agent. (NCI04); UMLS Semantic Type: STY:T200")
    assert node[NODE_DESCRIPTION] == "An agent. (NCI04)"
