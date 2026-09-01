"""Tests for qualifying LOINC Part names by their part type.

A LOINC Part is one axis-value of a LOINC term. The same name appears under several axes -- "Creatinine" is
both a COMPONENT (LP14355-9, the analyte) and a DIVISORS (LP32035-5, the denominator of a ratio) -- so the
part type is appended to disambiguate them.
"""

import pytest

from kraken.harmonizers.loinc import (
    PART_TYPE_COLUMN,
    LoincHarmonizer,
    NodeSpec,
    numeric_description,
    part_name_qualifier,
    part_type_display,
)


def _row(part_type, part_name=""):
    return {PART_TYPE_COLUMN: part_type, "PartName": part_name}


# ------------------------------- display term -------------------------------


def test_component_needs_no_qualifier():
    """COMPONENT is LOINC's unmarked default and 85% of parts -- a component part simply IS the analyte."""
    assert part_type_display("COMPONENT") is None


def test_missing_part_type_needs_no_qualifier():
    """Only Part.csv has the column; the other LOINC files pass through here with nothing set."""
    for absent in ("", "   ", None):
        assert part_type_display(absent) is None


def test_display_terms_are_lowercased():
    assert part_type_display("SYSTEM") == "system"
    assert part_type_display("METHOD") == "method"
    assert part_type_display("SUPER SYSTEM") == "super system"


def test_plural_part_type_is_singularised():
    """DIVISORS is LOINC's only plural. Handled by an explicit override rather than stripping a trailing 's',
    which would turn CLASS into 'clas'."""
    assert part_type_display("DIVISORS") == "divisor"
    assert part_type_display("CLASS") == "class"


def test_dotted_part_types_keep_their_family_and_axis():
    """First segment names the family, last is the axis; middle segments are elaboration. The family is what
    makes an otherwise-empty axis name mean something -- 'Discharge instructions [kind]' is a non-sequitur."""
    assert part_type_display("Document.Kind") == "document kind"
    assert part_type_display("Document.Role") == "document role"
    assert part_type_display("Rad.Guidance for.Action") == "rad action"
    assert part_type_display("Rad.Anatomic Location.Laterality.Presence") == "rad presence"


def test_middle_segments_are_dropped_avoiding_repetition():
    """A naive join of every segment would give "modality modality subtype" and "view view type"."""
    assert part_type_display("Rad.Modality.Modality Subtype") == "rad modality subtype"
    assert part_type_display("Rad.View.View Type") == "rad view type"


def test_camel_case_part_types_are_split():
    assert part_type_display("Document.SubjectMatterDomain") == "document subject matter domain"
    assert part_type_display("Document.TypeOfService") == "document type of service"


# ------------------------------- applying it to a name -------------------------------


def test_the_issues_two_creatinine_parts_become_distinguishable():
    """LP14355-9 and LP32035-5 are both named "Creatinine" in LOINC; only the divisor gets qualified."""
    assert part_name_qualifier(_row("COMPONENT", "Creatinine")) is None
    assert part_name_qualifier(_row("DIVISORS", "Creatinine")) == "divisor"


def test_names_that_already_say_their_type_are_left_alone():
    """LOINC's GENE parts are all named "<symbol> gene", and many SYSTEM parts contain "system" -- so this is
    idempotent and re-running can't double-qualify."""
    assert part_name_qualifier(_row("GENE", "AARS2 gene")) is None
    assert part_name_qualifier(_row("SYSTEM", "Respiratory system.airway+Inhl Gas")) is None
    assert part_name_qualifier(_row("METHOD", "Method X")) is None


def test_idempotence_check_is_case_insensitive():
    assert part_name_qualifier(_row("GENE", "AARS2 GENE")) is None


@pytest.mark.parametrize(
    "part_type,part_name,expected",
    [
        ("DIVISORS", "Cholesterol", "divisor"),
        ("CLASS", "PANEL.SURVEY.MFS", "class"),
        ("Rad.Anatomic Location.Imaging Focus", "Adrenal gland", "rad imaging focus"),
        ("Document.Kind", "Discharge instructions", "document kind"),
        ("Document.SubjectMatterDomain", "Forensic medicine", "document subject matter domain"),
        ("TIME", "30S", "time"),
    ],
)
def test_representative_parts_get_their_axis(part_type, part_name, expected):
    assert part_name_qualifier(_row(part_type, part_name)) == expected


# ------------------------------- qualifying a string -------------------------------


def test_qualify_appends_the_bracketed_term():
    assert LoincHarmonizer._qualify("Creatinine", "divisor") == "Creatinine [divisor]"


def test_qualify_is_idempotent():
    """Synonyms and names run through the same helper, and a re-run must not double-qualify."""
    once = LoincHarmonizer._qualify("Creatinine", "divisor")
    assert LoincHarmonizer._qualify(once, "divisor") == once
    assert LoincHarmonizer._qualify("Blood system", "system") == "Blood system"


# ------------------------------- descriptions -------------------------------


def test_definition_is_preferred_when_present():
    row = {"DefinitionDescription": "The disposition of an EMS unit.", "SURVEY_QUEST_TEXT": "Were you seen?"}
    assert numeric_description(row) == "The disposition of an EMS unit."


def test_survey_question_is_used_when_there_is_no_definition():
    """Labelled rather than bare: the prefix says what kind of thing the code is, which is the part that
    isn't already in the name."""
    row = {"DefinitionDescription": "", "SURVEY_QUEST_TEXT": "What is your race?"}
    assert numeric_description(row) == "Survey question: What is your race?"


def test_no_description_when_loinc_supplies_neither():
    """78% of numeric codes -- LOINC simply doesn't define most of them."""
    assert numeric_description({"DefinitionDescription": "", "SURVEY_QUEST_TEXT": ""}) is None
    assert numeric_description({}) is None


def test_description_spec_reads_a_named_column():
    """AnswerList.csv names its column directly rather than needing a callable."""
    spec = NodeSpec("AnswerStringId", "DisplayText", "biolink:InformationContentEntity", "Description")
    assert LoincHarmonizer._description({"Description": " No pain "}, spec) == "No pain"
    assert LoincHarmonizer._description({"Description": ""}, spec) is None


def test_specs_without_a_description_get_none():
    """Parts have no description column in LOINC at all."""
    spec = NodeSpec("PartNumber", "PartName", "biolink:NamedThing")
    assert LoincHarmonizer._description({"PartName": "Creatinine"}, spec) is None
