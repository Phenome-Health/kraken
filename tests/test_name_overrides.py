"""Tests for name overrides (issue #12).

Some vocabularies name a node after what it's *about* rather than what it *is* -- KEGG's Parkinson disease
pathway is called "Parkinson disease" -- which makes it indistinguishable from the real disease in a result
list. These cover the rule table itself and its integration with create_node, where the interaction with the
synonym auto-add matters.
"""

import pytest

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.name_overrides import (
    NAME_OVERRIDES,
    ORIGINAL_NAME_ATTRIBUTE,
    apply_name_override,
)
from kraken.utils.constants import NODE_ATTRIBUTES, NODE_NAME, NODE_SYNONYMS


class _StubBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


class _Harmonizer(BaseHarmonizer):
    source_infores = "infores:test"


@pytest.fixture
def harmonizer() -> BaseHarmonizer:
    """Allocated without __init__, which reaches the network for the Biolink model and biomapper2.

    normalize_curie is stubbed too: some prefixes under test (KEGG) route it through biomapper2, and curie
    normalization is a separate concern from naming."""
    instance = object.__new__(_Harmonizer)
    instance.biolink = _StubBiolink()
    instance.normalize_curie = lambda curie: curie
    instance.name_override_count = 0
    instance.multi_taxon_node_count = 0
    instance.multi_taxon_examples = []
    return instance


def _node(harmonizer, curie, categories, name, synonyms=None):
    return harmonizer.create_node(
        curie=curie, categories=categories, provided_by="infores:test", name=name, synonyms=synonyms
    )


# ------------------------------- the rule table -------------------------------


def test_the_issues_own_examples_are_fixed():
    """KEGG:05012 and LOINC:LA27533-1 are the nodes named in issue #12; Rhea is from its comment."""
    assert apply_name_override("KEGG:05012", ["biolink:Pathway"], "Parkinson disease") == "Parkinson disease pathway"
    assert (
        apply_name_override("LOINC:LA27533-1", ["biolink:InformationContentEntity"], "Parkinson's disease")
        == "Parkinson's disease [answer]"
    )
    assert apply_name_override("RHEA:31979", ["biolink:MolecularActivity"], "Creatininase") == "Creatininase reaction"


def test_overrides_are_idempotent():
    """Names that already say what they are must not be double-suffixed, and re-running must be safe."""
    for name in ("Signaling by WNT pathway", "Signaling by WNT PATHWAY"):
        assert apply_name_override("REACT:R-HSA-1", ["biolink:Pathway"], name) == name
    once = apply_name_override("KEGG:05012", ["biolink:Pathway"], "Parkinson disease")
    assert apply_name_override("KEGG:05012", ["biolink:Pathway"], once) == once


def test_rules_require_both_the_id_prefix_and_the_category():
    """The category alone is too blunt (GO MolecularActivity names mostly self-identify already) and the
    prefix alone too broad (LOINC's LP parts and LL lists aren't answers)."""
    # right prefix, wrong category
    assert apply_name_override("KEGG:05012", ["biolink:Disease"], "Parkinson disease") == "Parkinson disease"
    # right category, wrong prefix -- the actual disease node must never be touched
    assert apply_name_override("UMLS:C0030567", ["biolink:Disease"], "Parkinson's disease") == "Parkinson's disease"
    # LOINC parts share the prefix and category of LOINC answers but are keyed on "LOINC:LA"
    assert apply_name_override("LOINC:LP12345-6", ["biolink:InformationContentEntity"], "Sodium") == "Sodium"


def test_rhea_equations_are_left_alone_but_enzyme_names_are_not():
    """98% of Rhea names are chemical equations that already read as reactions; only the 9 named after their
    enzyme (the issue comment's "creatininase") need the suffix."""
    equation = "trimethylamine + NADPH + O2 = trimethylamine N-oxide"
    assert apply_name_override("RHEA:31979", ["biolink:MolecularActivity"], equation) == equation
    assert apply_name_override("RHEA:1", ["biolink:MolecularActivity"], "creatininase") == "creatininase reaction"
    assert apply_name_override("RHEA:2", ["biolink:MolecularActivity"], "phytol kinase") == "phytol kinase reaction"


def test_every_rule_appends_only():
    """Append-only is what lets the bare name be dropped from synonyms: the original stays a substring, so
    text search still finds it. A rule that rewrote the name would break that."""
    for (id_prefix, category), override in NAME_OVERRIDES.items():
        result = apply_name_override(f"{id_prefix}X", [category], "Some Name")
        assert result.startswith("Some Name"), f"{id_prefix}/{category} is not append-only"
        assert result != "Some Name", f"{id_prefix}/{category} did not apply"


# ------------------------------- integration with create_node -------------------------------


def test_synonym_gets_the_renamed_form_not_the_bare_name(harmonizer):
    """The auto-added synonym must be the renamed name. Keeping the bare form would rank the answer code
    alongside the real disease -- exactly the confusion the override removes."""
    node = _node(harmonizer, "LOINC:LA27533-1", ["biolink:InformationContentEntity"], "Parkinson's disease")
    assert node[NODE_NAME] == "Parkinson's disease [answer]"
    assert node[NODE_SYNONYMS] == ["Parkinson's disease [answer]"]


def test_original_name_is_kept_in_attributes(harmonizer):
    node = _node(harmonizer, "KEGG:05012", ["biolink:Pathway"], "Parkinson disease")
    assert node[NODE_ATTRIBUTES]["infores:test"][ORIGINAL_NAME_ATTRIBUTE] == "Parkinson disease"
    assert harmonizer.name_override_count == 1


def test_an_explicitly_passed_synonym_is_left_alone(harmonizer):
    """Only the auto-added synonym is affected. A harmonizer that deliberately supplies the raw name knows its
    source better than this rule does."""
    node = _node(harmonizer, "KEGG:05012", ["biolink:Pathway"], "Parkinson disease", synonyms=["Parkinson disease"])
    assert node[NODE_NAME] == "Parkinson disease pathway"
    assert set(node[NODE_SYNONYMS]) == {"Parkinson disease", "Parkinson disease pathway"}


def test_untouched_nodes_keep_their_name_and_gain_no_attribute(harmonizer):
    node = _node(harmonizer, "UMLS:C0030567", ["biolink:Disease"], "Parkinson's disease")
    assert node[NODE_NAME] == "Parkinson's disease"
    assert NODE_ATTRIBUTES not in node
    assert harmonizer.name_override_count == 0
