"""icd10cm: CDC's fixed-width file is read right, the code tree is real, and Mondo's rows can't merge what isn't one."""

import json
from types import SimpleNamespace

import jsonlines
import pytest

from kraken.entity_resolution.build import resolve_entities
from kraken.entity_resolution.families import BranchFamilies
from kraken.harmonizers.icd10cm import (
    CATEGORY_BY_CHAPTER,
    DISEASE_CATEGORY,
    ICD10CMHarmonizer,
    dotted,
    harmonized_titles,
    parent_code,
    read_order_file,
)
from tests.helpers import stub_normalization

SOURCE_ID = "infores:icd10cm"
SAME_AS = "biolink:same_as"
CLOSE_MATCH = "biolink:close_match"
SUBCLASS_OF = "biolink:subclass_of"


def _line(number: int, code: str, flag: int, name: str) -> str:
    return f"{number:05d} {code:<7} {flag} {name[:60]:<60} {name}\r\n"


ORDER_ROWS = [
    ("A08", 0, "Viral and other specified intestinal infections"),
    ("A083", 0, "Other viral enteritis"),
    ("A0839", 1, "Other viral enteritis"),
    ("D70", 0, "Neutropenia"),
    ("D82", 0, "Immunodeficiency associated with other major defects"),
    ("D828", 1, "Immunodeficiency associated with other specified major defects"),
    ("E11", 0, "Type 2 diabetes mellitus"),
    ("E116", 0, "Type 2 diabetes mellitus with other specified complications"),
    ("E1165", 1, "Type 2 diabetes mellitus with hyperglycemia"),
    ("E13", 0, "Other specified diabetes mellitus"),
    ("E78", 0, "Disorders of lipoprotein metabolism and other lipidemias"),
    ("E785", 1, "Hyperlipidemia, unspecified"),
    ("F53", 0, "Mental and behavioral disorders associated with the puerperium, not elsewhere classified"),
    ("F531", 1, "Puerperal psychosis"),
    ("I10", 1, "Essential (primary) hypertension"),
    ("I15", 0, "Secondary hypertension"),
    ("I21", 0, "Acute myocardial infarction"),
    ("I219", 1, "Acute myocardial infarction, unspecified"),
    ("I150", 1, "Renovascular hypertension"),
    ("J94", 0, "Other pleural conditions"),
    ("J942", 1, "Hemothorax"),
    ("K28", 0, "Gastrojejunal ulcer"),
    ("K280", 1, "Acute gastrojejunal ulcer with hemorrhage"),
    ("K76", 0, "Other diseases of liver"),
    ("K768", 0, "Other specified diseases of liver"),
    ("N18", 0, "Chronic kidney disease (CKD)"),
    ("N185", 1, "Chronic kidney disease, stage 5"),
    ("N186", 1, "End stage renal disease"),
    ("Q80", 0, "Congenital ichthyosis"),
    ("Q802", 1, "Lamellar ichthyosis"),
    ("R10", 0, "Abdominal and pelvic pain"),
    ("T36", 0, "Poisoning by, adverse effect of and underdosing of systemic antibiotics"),
    ("T360", 0, "Poisoning by, adverse effect of and underdosing of penicillins"),
    ("T360X1", 0, "Poisoning by penicillins, accidental (unintentional)"),
    ("T360X1A", 1, "Poisoning by penicillins, accidental (unintentional), initial encounter"),
    ("V00", 0, "Pedestrian conveyance accident"),
    ("Z15", 0, "Genetic susceptibility to disease"),
    ("Z79", 0, "Long term (current) drug therapy"),
]
ORDER_FILE = "".join(_line(i + 1, code, flag, name) for i, (code, flag, name) in enumerate(ORDER_ROWS))

MONDO_HEADER = "subject_id\tsubject_label\tpredicate_id\tobject_id\tobject_label\tmapping_justification\n"
EXACT = "skos:exactMatch"
MONDO_ROWS = [
    ("MONDO:0005148", "type 2 diabetes mellitus", EXACT, "ICD10CM:E11", "Type 2 diabetes mellitus"),
    ("MONDO:0001155", "gastrojejunal ulcer", EXACT, "ICD10CM:K28", "Gastrojejunal ulcer"),
    ("MONDO:0001155", "gastrojejunal ulcer", EXACT, "ICD10CM:K28.0", "Acute gastrojejunal ulcer with hemorrhage"),
    ("MONDO:0018623", "postpartum psychosis", EXACT, "ICD10CM:F53", "Mental and behavioral disorders"),
    ("MONDO:0018623", "postpartum psychosis", EXACT, "ICD10CM:F53.1", "Puerperal psychosis"),
    ("MONDO:0004375", "end stage renal failure", EXACT, "ICD10CM:N18.5", "Chronic kidney disease, stage 5"),
    ("MONDO:0004375", "end stage renal failure", EXACT, "ICD10CM:N18.6", "End stage renal disease"),
    ("MONDO:0005044", "hypertensive disorder", EXACT, "ICD10CM:I10-I16", ""),  # the real row's label is blank
    ("MONDO:0005044", "hypertensive disorder", EXACT, "ICD10CM:I15", "Secondary hypertension"),
    ("MONDO:7770003", "cystic fibrosis-related diabetes", EXACT, "ICD10CM:E13", "Other specified diabetes mellitus"),
    ("MONDO:0001609", "agranulocytosis", EXACT, "ICD10CM:D70", "Neutropenia"),
    ("MONDO:0001475", "neutropenia", EXACT, "DOID:1227", "neutropenia"),
    ("MONDO:0006783", "hemopneumothorax", EXACT, "ICD10CM:J94.2", "Hemothorax"),
    (
        "MONDO:0004781",
        "acute myocardial infarction",
        EXACT,
        "ICD10CM:I21.9",
        "Acute myocardial infarction, unspecified",
    ),
    ("MONDO:0020573", "inherited disease susceptibility", EXACT, "ICD10CM:Z15", "Genetic susceptibility to disease"),
    # none of these link anything: a range alone, an obsolete term, a broad row, another code system
    ("MONDO:0000916", "intestinal infectious disease", EXACT, "ICD10CM:A00-A09", "Intestinal infectious diseases"),
    ("MONDO:0001946", "obsolete hyperestrogenism", EXACT, "ICD10CM:E78", "Estrogen excess"),
    ("MONDO:0035511", "ricin poisoning", "skos:broadMatch", "ICD10CM:T36", "Toxic effect"),
    ("MONDO:0000367", "taeniasis", EXACT, "ICD10WHO:B68", "Taeniasis"),
]
MONDO_FILE = (
    "# curie_map:\n#   ICD10CM: https://example.org/\n"
    + MONDO_HEADER
    + "".join("\t".join((*row, "semapv:UnspecifiedMatching")) + "\n" for row in MONDO_ROWS)
)

# Real terms of the 2026-10-29 release, cut to the lines that matter
DOID_FILE = """format-version: 1.2
data-version: releases/2026-10-29/doid.obo

[Term]
id: DOID:10763
name: hypertension
xref: ICD10CM:I10

[Term]
id: DOID:10825
name: essential hypertension
xref: ICD10CM:I10

[Term]
id: DOID:1168
name: familial hyperlipidemia
xref: ICD10CM:E78.5

[Term]
id: DOID:9352
name: type 2 diabetes mellitus
xref: ICD10CM:E11

[Term]
id: DOID:1591
name: renovascular hypertension
xref: ICD10CM:I15.0

[Term]
id: DOID:0090113
name: RIDDLE syndrome
xref: ICD10CM:D82.8

[Term]
id: DOID:0060716
name: obsolete autosomal recessive congenital ichthyosis 7
xref: ICD10CM:Q80.2
is_obsolete: true

[Typedef]
id: has_symptom
name: has_symptom
"""


class _LeafBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


def _harmonizer() -> ICD10CMHarmonizer:
    """A harmonizer allocated without __init__ (which builds a Biolink toolkit and biomapper2 Normalizer)."""
    harmonizer = object.__new__(ICD10CMHarmonizer)
    harmonizer.source_infores = SOURCE_ID
    harmonizer.biolink = _LeafBiolink()
    harmonizer.name_override_count = 0
    harmonizer.multi_taxon_node_count = 0
    harmonizer.multi_taxon_examples = []
    return stub_normalization(harmonizer)


def _write_inputs(directory) -> None:
    directory.mkdir(exist_ok=True)
    (directory / "icd10cm-order-2027.txt").write_bytes(ORDER_FILE.encode())
    (directory / "mondo.sssom.2026-11-06.tsv").write_text(MONDO_FILE)
    (directory / "doid.2026-10-29.obo").write_text(DOID_FILE)


def _run(tmp_path) -> tuple[dict, dict]:
    """``id -> node`` and ``(subject, predicate, object) -> edge``."""
    inputs = tmp_path / "inputs"
    _write_inputs(inputs)
    nodes_path, edges_path = tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl"
    _harmonizer().harmonize(nodes_path, edges_path, input_file=inputs)
    nodes = {n["id"]: n for n in map(json.loads, nodes_path.read_text().splitlines())}
    edges = {
        (e["subject"], e["predicate"], e["object"]): e for e in map(json.loads, edges_path.read_text().splitlines())
    }
    return nodes, edges


def test_the_order_file_is_read_by_its_columns(tmp_path):
    # Two real lines, the second with a long name that differs from its short one
    path = tmp_path / "order.txt"
    path.write_bytes(
        b"00001 A00     0 Cholera                                                      Cholera\r\n"
        b"73123 T360X1A 1 Poisoning by penicillins, accidental (unintentional), init   "
        b"Poisoning by penicillins, accidental (unintentional), initial encounter\r\n"
    )
    entries = read_order_file(path)
    assert entries["A00"].billable is False and entries["A00"].name == "Cholera"
    assert entries["T360X1A"].billable is True
    assert entries["T360X1A"].name == "Poisoning by penicillins, accidental (unintentional), initial encounter"


def test_a_changed_layout_fails_the_build(tmp_path):
    path = tmp_path / "order.txt"
    path.write_text("00001 A00      0 Cholera   Cholera\n")
    with pytest.raises(ValueError, match="not an order file line"):
        read_order_file(path)


def test_the_dot_goes_after_the_third_character():
    assert dotted("E11") == "E11"
    assert dotted("E1165") == "E11.65"
    assert dotted("T360X1A") == "T36.0X1A"


def test_a_parent_is_the_longest_prefix_the_file_lists():
    codes = {"T36": None, "T360": None, "T360X1": None, "T360X1A": None}
    assert parent_code("T360X1A", codes) == "T360X1"
    assert parent_code("T360X1", codes) == "T360"  # there is no T360X: X only holds a place
    assert parent_code("T36", codes) is None
    with pytest.raises(ValueError, match="no parent"):
        parent_code("S72001", codes)


def test_codes_and_their_tree(tmp_path):
    nodes, edges = _run(tmp_path)
    assert nodes["ICD10CM:E11.65"]["attributes"][SOURCE_ID]["billable"] is True
    assert nodes["ICD10CM:E11"]["attributes"][SOURCE_ID]["billable"] is False
    assert ("ICD10CM:E11.65", SUBCLASS_OF, "ICD10CM:E11.6") in edges
    assert ("ICD10CM:T36.0X1", SUBCLASS_OF, "ICD10CM:T36.0") in edges
    assert not [e for e in edges if e[0] == "ICD10CM:E11" and e[1] == SUBCLASS_OF]


def test_mondo_exact_rows_are_same_as_with_mondo_as_subject(tmp_path):
    nodes, edges = _run(tmp_path)
    assert ("MONDO:0005148", SAME_AS, "ICD10CM:E11") in edges
    assert nodes["MONDO:0005148"]["categories"] == [DISEASE_CATEGORY]


def test_ranges_obsolete_terms_broad_rows_and_other_code_systems_link_nothing(tmp_path):
    nodes, edges = _run(tmp_path)
    assert not [n for n in nodes if "-" in n]
    for absent in ("MONDO:0000916", "MONDO:0001946", "MONDO:0035511", "MONDO:0000367"):
        assert absent not in nodes
    assert not [e for e in edges if e[0].startswith("MONDO:") and e[2] in ("ICD10CM:E78", "ICD10CM:T36")]


@pytest.mark.parametrize(
    "mondo_id, code",
    [
        ("MONDO:0001155", "ICD10CM:K28"),  # several codes for one term: the parent and child
        ("MONDO:0001155", "ICD10CM:K28.0"),
        ("MONDO:0018623", "ICD10CM:F53.1"),
        ("MONDO:0004375", "ICD10CM:N18.6"),  # two siblings
        ("MONDO:0005044", "ICD10CM:I15"),  # a code and a range
        ("MONDO:7770003", "ICD10CM:E13"),  # "Other specified diabetes mellitus"
        ("MONDO:0001609", "ICD10CM:D70"),  # "Neutropenia", while MONDO:0001475 is named neutropenia
        ("MONDO:0006783", "ICD10CM:J94.2"),  # reviewed: hemothorax is broader than hemopneumothorax
    ],
)
def test_a_row_mondo_is_not_clear_on_links_without_merging(tmp_path, mondo_id, code):
    nodes, edges = _run(tmp_path)
    assert (mondo_id, CLOSE_MATCH, code) in edges
    assert (mondo_id, SAME_AS, code) not in edges
    assert nodes[code]["name"]  # so a name match can still place it


def test_only_a_code_with_a_same_as_carries_no_name(tmp_path):
    nodes, _edges = _run(tmp_path)
    assert "name" not in nodes["ICD10CM:E11"]
    assert nodes["ICD10CM:E11"]["unmatched_name"] == "Type 2 diabetes mellitus"
    assert nodes["ICD10CM:E11"]["synonyms"] == ["Type 2 diabetes mellitus"]
    assert nodes["ICD10CM:E11.65"]["name"] == "Type 2 diabetes mellitus with hyperglycemia"


def test_chapters_decide_the_category_unless_mondo_says_disease(tmp_path):
    nodes, _edges = _run(tmp_path)
    assert nodes["ICD10CM:K28.0"]["categories"] == ["biolink:Disease"]
    assert nodes["ICD10CM:R10"]["categories"] == ["biolink:PhenotypicFeature"]
    assert nodes["ICD10CM:T36"]["categories"] == ["biolink:DiseaseOrPhenotypicFeature"]
    assert nodes["ICD10CM:V00"]["categories"] == ["biolink:Event"]
    assert nodes["ICD10CM:Z79"]["categories"] == ["biolink:ClinicalEntity"]
    assert nodes["ICD10CM:Z15"]["categories"] == ["biolink:Disease"]  # or the family guardrail would cut its same_as


def test_every_category_sits_in_one_family():
    """A bridge category such as ClinicalFinding is not inherited by entity resolution, so an id Babel doesn't know
    falls through to NamedThing, a wildcard that can merge on name with anything."""
    families = BranchFamilies.load()
    for category in {DISEASE_CATEGORY, *CATEGORY_BY_CHAPTER.values()}:
        assert len(families.families_for_category(category)) == 1, category


def test_disease_ontology_links_a_code_without_its_own_mondo_row(tmp_path):
    nodes, edges = _run(tmp_path)
    assert nodes["DOID:10825"]["name"] == "essential hypertension"  # so a term newer than our Babel isn't nameless
    for doid_id in ("DOID:10763", "DOID:10825"):  # I10 is listed by both terms, and both are linked
        assert edges[(doid_id, CLOSE_MATCH, "ICD10CM:I10")]["primary_knowledge_source"] == "infores:disease-ontology"
    assert ("DOID:1168", CLOSE_MATCH, "ICD10CM:E78.5") in edges
    assert ("DOID:1591", CLOSE_MATCH, "ICD10CM:I15.0") in edges  # its parent I15 has a Mondo row, it doesn't
    assert not [e for e in edges if e[0].startswith("DOID:") and e[1] != CLOSE_MATCH]


@pytest.mark.parametrize(
    "doid_id",
    [
        "DOID:9352",  # E11 has its own Mondo row
        "DOID:0090113",  # D82.8 "... associated with other specified major defects" is a residual group
        "DOID:0060716",  # an obsolete term
    ],
)
def test_disease_ontology_links_only_real_codes_mondo_leaves(tmp_path, doid_id):
    nodes, edges = _run(tmp_path)
    assert doid_id not in nodes
    assert not [e for e in edges if doid_id in e]


# --- through the real entity resolution ------------------------------------------------------------------------------


def _babel_node(curie: str, name: str) -> dict:
    return {
        "id": curie,
        "name": name,
        "categories": [DISEASE_CATEGORY],
        "provided_by": ["infores:sri-node-normalizer"],
        "equivalent_ids": [curie],
    }


def _babel_clique(hub: str, *others: str) -> list[dict]:
    return [
        {
            "subject": hub,
            "predicate": SAME_AS,
            "object": other,
            "primary_knowledge_source": "infores:sri-node-normalizer",
        }
        for other in others
    ]


def _resolve(tmp_path) -> dict[str, dict]:
    """``id -> its merged node`` from real ids: T2D's clique, UMLS's namesakes of CM codes, and a phenotype named like
    a code that Mondo calls exact to another disease."""
    harmonized = tmp_path / "harmonized"
    (harmonized / "babel").mkdir(parents=True)
    (harmonized / "icd10cm").mkdir()
    babel_nodes = [
        _babel_node("MONDO:0005148", "type 2 diabetes mellitus"),
        _babel_node("DOID:9352", "type 2 diabetes mellitus"),
        _babel_node("NCIT:C195403", "Type 2 Diabetes Mellitus with Hyperglycemia"),
        _babel_node("UMLS:C0494062", "Other viral enteritis"),
        _babel_node("MONDO:0001155", "gastrojejunal ulcer"),
        _babel_node("MONDO:0006783", "hemopneumothorax"),
        _babel_node("HP:0012151", "Hemothorax"),
        _babel_node("MONDO:0001475", "neutropenia"),
        _babel_node("MONDO:0001609", "agranulocytosis"),
        _babel_node("MONDO:0004781", "acute myocardial infarction"),
        _babel_node("DOID:9408", "acute myocardial infarction"),
        _babel_node("NCIT:C35204", "Acute Myocardial Infarction"),
        _babel_node("UMLS:C0155626", "Acute myocardial infarction"),
    ]
    with jsonlines.open(harmonized / "babel" / "nodes.jsonl", "w") as w:
        w.write_all(babel_nodes)
    with jsonlines.open(harmonized / "babel" / "edges.jsonl", "w") as w:
        w.write_all(_babel_clique("MONDO:0005148", "DOID:9352"))
        w.write_all(_babel_clique("MONDO:0004781", "DOID:9408", "NCIT:C35204", "UMLS:C0155626"))

    inputs = tmp_path / "inputs"
    _write_inputs(inputs)
    nodes_path, edges_path = harmonized / "icd10cm" / "nodes.jsonl", harmonized / "icd10cm" / "edges.jsonl"
    _harmonizer().harmonize(nodes_path, edges_path, input_file=inputs)

    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={
            "babel": (harmonized / "babel" / "nodes.jsonl", harmonized / "babel" / "edges.jsonl"),
            "icd10cm": (nodes_path, edges_path),
        },
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    resolve_entities(config, biolink=None)
    with jsonlines.open(config.integrated_nodes_path) as reader:
        return {member: node for node in reader for member in node["equivalent_ids"]}


def test_a_three_character_code_joins_its_disease_and_its_child_stays_out(tmp_path):
    by_id = _resolve(tmp_path)
    t2d = by_id["MONDO:0005148"]
    assert "ICD10CM:E11" in t2d["equivalent_ids"]
    assert t2d["name"] == "type 2 diabetes mellitus"
    assert "Type 2 diabetes mellitus" in t2d["synonyms"]
    assert "ICD10CM:E11.65" not in t2d["equivalent_ids"]
    # ...and where a source has E11.65's own concept, the code goes there
    assert by_id["ICD10CM:E11.65"] is by_id["NCIT:C195403"]


def test_a_parent_and_child_with_one_name_stay_split(tmp_path):
    by_id = _resolve(tmp_path)
    assert by_id["ICD10CM:A08.3"] is not by_id["ICD10CM:A08.39"]
    assert by_id["UMLS:C0494062"]["equivalent_ids"] == ["UMLS:C0494062"]  # the namesake of both: a tie, so neither


def test_a_term_with_two_codes_merges_only_the_one_named_like_it(tmp_path):
    by_id = _resolve(tmp_path)
    assert by_id["ICD10CM:K28"] is by_id["MONDO:0001155"]
    assert by_id["ICD10CM:K28.0"] is not by_id["MONDO:0001155"]


def test_a_code_named_like_another_term_goes_to_that_term(tmp_path):
    by_id = _resolve(tmp_path)
    assert by_id["ICD10CM:D70"] is by_id["MONDO:0001475"]
    assert by_id["ICD10CM:D70"] is not by_id["MONDO:0001609"]


def test_a_code_broader_than_its_mondo_term_stays_out_of_it(tmp_path):
    """Mondo calls J94.2 "Hemothorax" exact to hemopneumothorax, a row reviewed as wrong. Its name places it."""
    by_id = _resolve(tmp_path)
    assert by_id["ICD10CM:J94.2"] is not by_id["MONDO:0006783"]
    assert by_id["ICD10CM:J94.2"] is by_id["HP:0012151"]


def test_a_code_its_disease_has_no_room_for_is_named_by_its_title(tmp_path):
    """Mondo calls I21.9 exact to acute myocardial infarction, so the code carries no name to match on. But I21, named
    like the disease, is in it, and a node holds one ICD-10-CM code. Eight such codes had no name at all in a build."""
    by_id = _resolve(tmp_path)
    assert by_id["ICD10CM:I21"] is by_id["MONDO:0004781"]
    left_out = by_id["ICD10CM:I21.9"]
    assert left_out["equivalent_ids"] == ["ICD10CM:I21.9"]
    assert left_out["name"] == "Acute myocardial infarction, unspecified"


def test_titles_read_back_from_the_harmonized_nodes_are_the_order_files(tmp_path):
    """The icd10 source reads them there, whether or not the name is one to match on."""
    inputs = tmp_path / "inputs"
    _write_inputs(inputs)
    nodes_path = tmp_path / "nodes.jsonl"
    _harmonizer().harmonize(nodes_path, tmp_path / "edges.jsonl", input_file=inputs)
    assert harmonized_titles(nodes_path) == {code: name for code, _billable, name in ORDER_ROWS}
