"""icd10: WHO's code file is read right, the code tree is real, a code Babel already places is left alone, and a code
is its ICD-10-CM twin only when both have one title."""

import json
import shutil
from types import SimpleNamespace

import jsonlines
import pytest

from kraken.entity_resolution.build import resolve_entities
from kraken.harmonizers.icd10 import ICD10Harmonizer, code_names, output_is_current, read_codes_file, title_key
from tests.helpers import stub_normalization

SOURCE_ID = "who-icd10"
SUBCLASS_OF = "biolink:subclass_of"
SAME_AS = "biolink:same_as"
DISEASE = "biolink:Disease"

# Real lines of the 2019 edition, without their last five fields (references to WHO's tabulation lists)
CODE_ROWS = [
    "3;N;X;04;E10;E11.-;E11;E11;Type 2 diabetes mellitus;Type 2 diabetes mellitus;;",
    "4;T;S;04;E10;E11.0;E11.0;E110;Type 2 diabetes mellitus: With coma;Type 2 diabetes mellitus;With coma;",
    "3;T;X;09;I10;I10;I10;I10;Essential (primary) hypertension;Essential (primary) hypertension;;",
    "3;N;X;09;I20;I21.-;I21;I21;Acute myocardial infarction;Acute myocardial infarction;;",
    "4;T;X;09;I20;I21.9;I21.9;I219;Acute myocardial infarction, unspecified;Acute myocardial infarction;Acute "
    "myocardial infarction, unspecified;",
    "3;N;X;01;B15;B18.-;B18;B18;Chronic viral hepatitis;Chronic viral hepatitis;;",
    "4;N;X;01;B15;B18.0;B18.0;B180;Chronic viral hepatitis B with delta-agent;Chronic viral hepatitis;Chronic viral "
    "hepatitis B with delta-agent;",
    "5;T;S;01;B15;B18.00;B18.00;B1800;Chronic viral hepatitis B with delta-agent: immune-tolerant phase;Chronic viral "
    "hepatitis;Chronic viral hepatitis B with delta-agent;immune-tolerant phase",
    "3;N;X;03;D60;D63.-*;D63;D63;Anaemia in chronic diseases classified elsewhere;Anaemia in chronic diseases "
    "classified elsewhere;;",
    "3;N;X;03;D60;D64.-;D64;D64;Other anaemias;Other anaemias;;",
    "4;T;X;03;D60;D64.9;D64.9;D649;Anaemia, unspecified;Other anaemias;Anaemia, unspecified;",
    "3;N;X;06;G50;G50.-;G50;G50;Disorders of trigeminal nerve;Disorders of trigeminal nerve;;",
    "4;T;X;06;G50;G50.1;G50.1;G501;Atypical facial pain;Disorders of trigeminal nerve;Atypical facial pain;",
    "3;N;X;06;G60;G61.-;G61;G61;Inflammatory polyneuropathy;Inflammatory polyneuropathy;;",
    "4;T;X;06;G60;G61.0;G61.0;G610;Guillain-Barré syndrome;Inflammatory polyneuropathy;Guillain-Barré syndrome;",
    "3;N;X;17;Q90;Q93.-;Q93;Q93;Monosomies and deletions from the autosomes, not elsewhere classified;Monosomies and "
    "deletions from the autosomes, not elsewhere classified;;",
    "4;T;X;17;Q90;Q93.5;Q93.5;Q935;Other deletions of part of a chromosome;Monosomies and deletions from the "
    "autosomes, not elsewhere classified;Other deletions of part of a chromosome;",
    "3;N;X;18;R10;R10.-;R10;R10;Abdominal and pelvic pain;Abdominal and pelvic pain;;",
    "3;N;X;20;X85;Y06.-;Y06;Y06;Neglect and abandonment;Neglect and abandonment;;",
    "4;T;X;20;X85;Y06.0;Y06.0;Y060;By spouse or partner;Neglect and abandonment;By spouse or partner;",
    "3;N;X;20;X85;Y07.-;Y07;Y07;Other maltreatment;Other maltreatment;;",
    "4;T;X;20;X85;Y07.0;Y07.0;Y070;By spouse or partner;Other maltreatment;By spouse or partner;",
    "3;N;X;20;Y40;Y44.-;Y44;Y44;Agents primarily affecting blood constituents;Agents primarily affecting blood "
    "constituents;;",
    "4;T;X;20;Y40;Y44.2;Y44.2;Y442;Anticoagulants;Agents primarily affecting blood constituents;Anticoagulants;",
    "3;N;X;20;Y83;Y84.-;Y84;Y84;Other medical procedures as the cause of abnormal reaction of the patient, or of "
    "later complication, without mention of misadventure at the time of the procedure;Other medical procedures as "
    "the cause of abnormal reaction of the patient, or of later complication, without mention of misadventure at "
    "the time of the procedure;;",
    "4;T;X;20;Y83;Y84.0;Y84.0;Y840;Cardiac catheterization;Other medical procedures as the cause of abnormal "
    "reaction of the patient, or of later complication, without mention of misadventure at the time of the "
    "procedure;Cardiac catheterization;",
]
CODES_FILE = "".join(f"{row};001;4-001;3-001;2-001;1-001\r\n" for row in CODE_ROWS)

CM = "infores:icd10cm"
# The ICD-10-CM codes of those codes with CDC's titles, as the icd10cm source harmonizes them: the title is a code's
# name, one that is not matched on where the code merges through Mondo
ICD10CM_TITLES = {
    "B18": "Chronic viral hepatitis",
    "B18.0": "Chronic viral hepatitis B with delta-agent",
    "D64": "Other anemias",
    "D64.9": "Anemia, unspecified",
    "E11": "Type 2 diabetes mellitus",
    "E11.0": "Type 2 diabetes mellitus with hyperosmolarity",
    "G50": "Disorders of trigeminal nerve",
    "G50.1": "Atypical facial pain",
    "G61": "Inflammatory polyneuropathy",
    "G61.0": "Guillain-Barre syndrome",
    "I10": "Essential (primary) hypertension",
    "I21": "Acute myocardial infarction",
    "I21.9": "Acute myocardial infarction, unspecified",
    "Q93": "Monosomies and deletions from the autosomes, not elsewhere classified",
    "Q93.5": "Other deletions of part of a chromosome",
    "R10": "Abdominal and pelvic pain",
    "Y07": "Perpetrator of assault, maltreatment and neglect",
    "Y07.0": "Spouse or partner, perpetrator of maltreatment and neglect",
}
MONDO_TERM_OF_CODE = {"E11": "MONDO:0005148", "G50.1": "MONDO:0018362", "I21.9": "MONDO:0004781"}


def _icd10cm_nodes(titles: dict[str, str] = ICD10CM_TITLES) -> list[dict]:
    return [
        {
            "id": f"ICD10CM:{code}",
            ("unmatched_name" if code in MONDO_TERM_OF_CODE else "name"): title,
            "categories": [DISEASE],
            "provided_by": [CM],
            "equivalent_ids": [f"ICD10CM:{code}"],
            "synonyms": [title],
        }
        for code, title in titles.items()
    ]


BABEL = "infores:sri-node-normalizer"
# Babel has two of WHO's codes in a disease's clique, as cross-references: E11 on type 2 diabetes, which it is, and
# Q93.5 "Other deletions of part of a chromosome" on a syndrome it is far broader than
BABEL_EDGES = [
    {"subject": "MONDO:0005148", "predicate": SAME_AS, "object": "ICD10:E11", "primary_knowledge_source": BABEL},
    {"subject": "MONDO:0012315", "predicate": SAME_AS, "object": "ICD10:Q93.5", "primary_knowledge_source": BABEL},
]


class _LeafBiolink:
    """NamedThing, every category's ancestor, is no leaf beside another."""

    def filter_to_leaf_categories(self, categories):
        return [category for category in categories if category != "biolink:NamedThing"] or list(categories)


def _harmonizer(babel_edges, icd10cm_nodes) -> ICD10Harmonizer:
    """A harmonizer allocated without __init__ (which builds a Biolink toolkit and biomapper2 Normalizer)."""
    harmonizer = object.__new__(ICD10Harmonizer)
    harmonizer.source_infores = SOURCE_ID
    harmonizer.babel_edges = babel_edges
    harmonizer.icd10cm_nodes = icd10cm_nodes
    harmonizer.biolink = _LeafBiolink()
    harmonizer.name_override_count = 0
    harmonizer.multi_taxon_node_count = 0
    harmonizer.multi_taxon_examples = []
    return stub_normalization(harmonizer)


def _write(path, records) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with jsonlines.open(path, "w") as writer:
        writer.write_all(records)


def _harmonize(tmp_path, output_dir=None, babel_edges=None, codes_file: str = CODES_FILE) -> tuple:
    """Beside the icd10cm source's harmonized nodes, in ``tmp_path / "icd10cm"``."""
    inputs, icd10cm_nodes = tmp_path / "inputs", tmp_path / "icd10cm" / "nodes.jsonl"
    inputs.mkdir(exist_ok=True)
    (inputs / "icd102019syst_codes.txt").write_bytes(codes_file.encode("utf-8"))
    _write(icd10cm_nodes, _icd10cm_nodes())
    output_dir = output_dir or tmp_path
    nodes_path, edges_path = output_dir / "nodes.jsonl", output_dir / "edges.jsonl"
    _harmonizer(babel_edges, icd10cm_nodes).harmonize(nodes_path, edges_path, input_file=inputs)
    return nodes_path, edges_path


def _codes_file(tmp_path, text: str = CODES_FILE):
    path = tmp_path / "codes.txt"
    path.write_bytes(text.encode("utf-8"))
    return path


def test_whos_code_file_is_read_by_its_fields(tmp_path):
    entries = read_codes_file(_codes_file(tmp_path))
    assert entries["E110"].title == "Type 2 diabetes mellitus: With coma" and entries["E110"].leaf
    assert not entries["E11"].leaf
    assert "D63" in entries  # the code without its asterisk
    assert entries["G610"].title == "Guillain-Barré syndrome"  # the file is UTF-8
    assert entries["Y060"].category_title == "Neglect and abandonment"
    assert len(entries) == len(CODE_ROWS)


@pytest.mark.parametrize(
    ("text", "error"),
    [
        (CODES_FILE.replace(";001;4-001", ";4-001", 1), "not a code file line"),
        (CODES_FILE.replace("3;N;X;04;E10;E11.-", "4;N;X;04;E10;E11.-"), "not a code file line"),
        (CODES_FILE + CODES_FILE.splitlines(keepends=True)[0], "listed twice"),
        ("", "lists no code"),
    ],
    ids=["a field less", "the code is not as long as its level", "a code twice", "an empty file"],
)
def test_a_broken_file_fails_the_build(tmp_path, text, error):
    with pytest.raises(ValueError, match=error):
        read_codes_file(_codes_file(tmp_path, text))


def test_a_title_several_codes_share_goes_under_its_category(tmp_path):
    """WHO prints a subdivision under its category. Named "By spouse or partner" alone, neglect by a partner and
    maltreatment by a partner would be one node."""
    names = code_names(read_codes_file(_codes_file(tmp_path)))
    assert names["Y060"] == "Neglect and abandonment: By spouse or partner"
    assert names["Y070"] == "Other maltreatment: By spouse or partner"
    assert names["I10"] == "Essential (primary) hypertension"


def test_titles_compare_without_case_accents_and_british_spelling():
    assert title_key("Anaemia, unspecified") == title_key("Anemia, unspecified")
    assert title_key("Guillain-Barré syndrome") == title_key("Guillain-Barre syndrome")
    assert title_key("Labour and delivery complicated by cord complication") == title_key(
        "Labor and delivery complicated by cord complication"
    )
    assert title_key("Type 2 diabetes mellitus: With coma") != title_key(
        "Type 2 diabetes mellitus with hyperosmolarity"
    )


def _run(tmp_path) -> tuple[dict, dict]:
    """``id -> node`` and ``(subject, predicate, object) -> edge``, beside a Babel that places E11 and Q93.5."""
    babel_edges = tmp_path / "babel_edges.jsonl"
    _write(babel_edges, BABEL_EDGES)
    nodes_path, edges_path = _harmonize(tmp_path, babel_edges=babel_edges)
    nodes = {n["id"]: n for n in map(json.loads, nodes_path.read_text().splitlines())}
    edges = {
        (e["subject"], e["predicate"], e["object"]): e for e in map(json.loads, edges_path.read_text().splitlines())
    }
    return nodes, edges


def test_codes_and_their_tree(tmp_path):
    nodes, edges = _run(tmp_path)
    assert nodes["ICD10:E11.0"]["name"] == "Type 2 diabetes mellitus: With coma"
    assert nodes["ICD10:E11.0"]["attributes"][SOURCE_ID]["leaf"] is True
    assert nodes["ICD10:B18.0"]["attributes"][SOURCE_ID]["leaf"] is False
    assert nodes["ICD10:R10"]["categories"] == ["biolink:PhenotypicFeature"]  # ICD-10-CM's category for the chapter
    assert {key for key in edges if key[1] == SUBCLASS_OF} == {
        ("ICD10:B18.0", SUBCLASS_OF, "ICD10:B18"),
        ("ICD10:B18.00", SUBCLASS_OF, "ICD10:B18.0"),
        ("ICD10:D64.9", SUBCLASS_OF, "ICD10:D64"),
        ("ICD10:G50.1", SUBCLASS_OF, "ICD10:G50"),
        ("ICD10:G61.0", SUBCLASS_OF, "ICD10:G61"),
        ("ICD10:I21.9", SUBCLASS_OF, "ICD10:I21"),
        ("ICD10:Y06.0", SUBCLASS_OF, "ICD10:Y06"),
        ("ICD10:Y07.0", SUBCLASS_OF, "ICD10:Y07"),
        ("ICD10:Y44.2", SUBCLASS_OF, "ICD10:Y44"),
        ("ICD10:Y84.0", SUBCLASS_OF, "ICD10:Y84"),
    }


def test_a_code_babel_places_is_left_as_it_is(tmp_path):
    """Not a node, so no name goes to the disease Babel has it in, and no edge from it or to it hangs WHO's tree on
    that disease: neither Q93.5 under its category nor E11.0 under E11."""
    nodes, edges = _run(tmp_path)
    assert "ICD10:Q93.5" not in nodes and "ICD10:E11" not in nodes
    assert not [key for key in edges if {"ICD10:Q93.5", "ICD10:E11"} & set(key)]
    assert "ICD10:E11.0" in nodes


def test_an_icd10cm_twin_is_a_stub(tmp_path):
    """An edge needs both its nodes, and this one says nothing of the code the icd10cm source describes."""
    nodes, _edges = _run(tmp_path)
    stub = nodes["ICD10CM:I10"]
    assert stub.get("name") is None and stub["categories"] == ["biolink:NamedThing"]
    assert "attributes" not in stub and "synonyms" not in stub


def test_what_it_reads_has_to_be_harmonized_first(tmp_path):
    with pytest.raises(FileNotFoundError, match="Babel must be harmonized before"):
        _harmonize(tmp_path, babel_edges=tmp_path / "missing" / "edges.jsonl")
    inputs = tmp_path / "inputs"
    with pytest.raises(FileNotFoundError, match="icd10cm must be harmonized before"):
        _harmonizer(None, tmp_path / "missing" / "nodes.jsonl").harmonize(
            tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl", input_file=inputs
        )


def test_a_code_is_its_icd10cm_twin_only_with_the_same_title(tmp_path):
    nodes, edges = _run(tmp_path)
    twins = {key[0] for key in edges if key[1] == SAME_AS}
    assert twins == {
        f"ICD10:{code}"
        for code in (
            *("B18", "B18.0", "I10", "I21", "I21.9", "R10", "D64", "D64.9", "G50", "G50.1", "G61", "G61.0", "Q93"),
        )
    }
    # Its name is not one to match on, as for an ICD-10-CM code that merges through Mondo
    assert "name" not in nodes["ICD10:D64.9"] and nodes["ICD10:D64.9"]["unmatched_name"] == "Anaemia, unspecified"
    assert nodes["ICD10:D64.9"]["synonyms"] == ["Anaemia, unspecified"]
    assert nodes["ICD10:E11.0"]["name"] == "Type 2 diabetes mellitus: With coma"  # no twin: its name places it
    twin = edges[("ICD10:I10", SAME_AS, "ICD10CM:I10")]
    assert twin["primary_knowledge_source"] == "kraken", "our inference, not WHO's statement or CDC's"
    assert "aggregator_knowledge_source" not in twin, "and nothing we retrieved"
    assert twin["supporting_data_sources"] == [SOURCE_ID, "infores:icd10cm"]


def test_a_code_titled_by_the_drug_or_the_procedure_alone_is_not_matched_on(tmp_path):
    """Y44.2 is harm from anticoagulants in therapeutic use, and Y84.0 an abnormal reaction to a cardiac
    catheterization. Matched on its title, each would be the drug class or the procedure."""
    nodes, _edges = _run(tmp_path)
    titles = {
        "Y44": "Agents primarily affecting blood constituents",
        "Y44.2": "Anticoagulants",
        "Y84.0": "Cardiac catheterization",
    }
    for code, title in titles.items():
        assert "name" not in nodes[f"ICD10:{code}"] and nodes[f"ICD10:{code}"]["unmatched_name"] == title
    assert nodes["ICD10:Y84"]["name"].startswith("Other medical procedures as the cause of abnormal reaction")


def test_the_output_remembers_what_it_was_made_from(tmp_path):
    """A run often harmonizes one source. Output made beside one Babel must not be integrated beside another, or
    beside none: the codes that Babel placed are not in it, and the ones it did not place are named. Nor beside
    another ICD-10-CM list: its twins are that list's."""
    babel_edges, icd10cm_nodes = tmp_path / "babel_edges.jsonl", tmp_path / "icd10cm" / "nodes.jsonl"
    _write(babel_edges, BABEL_EDGES)
    nodes_path, _edges_path = _harmonize(tmp_path, babel_edges=babel_edges)
    assert output_is_current(nodes_path, babel_edges, icd10cm_nodes)
    copied = shutil.copy(babel_edges, tmp_path / "copied_edges.jsonl")
    assert output_is_current(nodes_path, copied, icd10cm_nodes), "the same Babel output, moved to another machine"

    assert not output_is_current(nodes_path, None, icd10cm_nodes), "a build without Babel"
    assert not output_is_current(nodes_path, babel_edges, None), "a build without ICD-10-CM"
    _write(copied, BABEL_EDGES[:1])
    assert not output_is_current(nodes_path, copied, icd10cm_nodes), "a Babel that does not place Q93.5"
    retitled = tmp_path / "retitled_nodes.jsonl"
    _write(retitled, _icd10cm_nodes({**ICD10CM_TITLES, "E11.0": "Type 2 diabetes mellitus with coma"}))
    assert not output_is_current(nodes_path, babel_edges, retitled), "a list in which E11.0 is WHO's E11.0"

    with pytest.raises(ValueError, match="lists no code"):
        _harmonize(tmp_path, babel_edges=babel_edges, codes_file="")
    assert not output_is_current(nodes_path, babel_edges, icd10cm_nodes), "a run that died on the way"


def _node(curie: str, name: str | None, source: str, category: str = DISEASE) -> dict:
    return {"id": curie, "name": name, "categories": [category], "provided_by": [source], "equivalent_ids": [curie]}


def _resolve(tmp_path, twin_statements: bool = True) -> dict[str, dict]:
    """``id -> its node`` after entity resolution, beside Babel's cliques, ICD-10-CM's codes as its harmonizer writes
    them, and a kg2 that lists ICD10:B18 on a Mondo term, types ICD10CM:I10 in two families, has a name of its own
    for ICD10:I21.9 and gives ICD10:G61.0 a species, which no ICD code has."""
    sources = ("babel", "icd10", "icd10cm", "kg2")
    paths = {s: (tmp_path / s / "nodes.jsonl", tmp_path / s / "edges.jsonl") for s in sources}
    babel_nodes = [
        _node("MONDO:0005148", "type 2 diabetes mellitus", BABEL),
        _node("ICD10:E11", None, BABEL),
        _node("MONDO:0012315", "distal 10q deletion syndrome", BABEL),
        _node("ICD10:Q93.5", None, BABEL),
        _node("MONDO:0018362", "persistent idiopathic facial pain", BABEL),
        _node("UMLS:C0154729", "Atypical facial pain", BABEL),
        _node("MONDO:0004781", "acute myocardial infarction", BABEL),
        _node("DOID:9408", "acute myocardial infarction", BABEL),
        _node("NCIT:C35204", "Acute Myocardial Infarction", BABEL),
        _node("UMLS:C0155626", "Acute myocardial infarction", BABEL),
        _node("UMLS:C0400918", "Chronic viral hepatitis B with hepatitis D", BABEL),
        _node("NCIT:C194012", "Chronic Viral Hepatitis B with Delta-Agent", BABEL),
        _node("SNOMEDCT:235869004", None, BABEL),
    ]
    cliques = [
        {"subject": hub, "predicate": SAME_AS, "object": other, "primary_knowledge_source": BABEL}
        for hub, others in (
            ("MONDO:0004781", ("DOID:9408", "NCIT:C35204", "UMLS:C0155626")),
            ("UMLS:C0400918", ("NCIT:C194012", "SNOMEDCT:235869004")),
        )
        for other in others
    ]
    mondo_rows = [
        {
            "subject": term,
            "predicate": SAME_AS,
            "object": f"ICD10CM:{code}",
            "primary_knowledge_source": "infores:mondo",
        }
        for code, term in MONDO_TERM_OF_CODE.items()
    ]
    animal_hepatitis = {
        **_node("MONDO:0025085", "hepatitis, viral, animal", "infores:rtx-kg2"),
        "equivalent_ids": ["MONDO:0025085", "ICD10:B18"],
    }
    _write(paths["babel"][0], babel_nodes)
    _write(paths["babel"][1], [*BABEL_EDGES, *cliques])
    _write(paths["icd10cm"][1], mondo_rows)
    hypertension = {**_node("ICD10CM:I10", None, "infores:rtx-kg2"), "categories": [DISEASE, "biolink:ClinicalEntity"]}
    infarction_nos = {
        **_node("ICD10:I21.9", "Acute myocardial infarction NOS", "infores:rtx-kg2"),
        "description": "UMLS Semantic Type: STY:T047",
    }
    of_a_mouse = {**_node("ICD10:G61.0", None, "infores:rtx-kg2"), "taxon": "NCBITaxon:10090"}
    _write(paths["kg2"][0], [animal_hepatitis, hypertension, infarction_nos, of_a_mouse])
    paths["icd10"][0].parent.mkdir(parents=True)
    nodes_path, edges_path = _harmonize(tmp_path, paths["icd10"][0].parent, babel_edges=paths["babel"][1])
    if not twin_statements:  # nor the stub that is a statement's other end
        _write(edges_path, [edge for edge in jsonlines.open(edges_path) if edge["predicate"] != SAME_AS])
        _write(nodes_path, [node for node in jsonlines.open(nodes_path) if node["id"].startswith("ICD10:")])

    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved=paths, integrated_dir=integrated, integrated_nodes_path=integrated / "nodes.jsonl"
    )
    resolve_entities(config, biolink=_LeafBiolink())
    with jsonlines.open(config.integrated_nodes_path) as reader:
        return {member: node for node in reader for member in node["equivalent_ids"]}


def test_codes_through_entity_resolution(tmp_path):
    by_id = _resolve(tmp_path)

    hypertension = by_id["ICD10:I10"]
    assert hypertension is by_id["ICD10CM:I10"]
    assert hypertension["id"] == "ICD10CM:I10", "a WHO code that joins a node must not rename it"
    assert by_id["ICD10:D64.9"] is by_id["ICD10CM:D64.9"], "anaemia is anemia"
    assert by_id["ICD10:E11.0"] is not by_id["ICD10CM:E11.0"], "with coma is not with hyperosmolarity"
    assert by_id["ICD10:Y06.0"] is not by_id["ICD10:Y07.0"], "neglect by a partner is not maltreatment by a partner"
    # B18.0 sits in a concept of three ids on one name match, NCIT's. Its WHO twin joins it there and must not take it
    # out: weighed as evidence, their same_as merged the two first, and one name match does not hold a group of two
    hepatitis = by_id["NCIT:C194012"]
    assert by_id["ICD10CM:B18.0"] is hepatitis and by_id["ICD10:B18.0"] is hepatitis
    assert {"UMLS:C0400918", "SNOMEDCT:235869004"} <= set(hepatitis["equivalent_ids"])
    # B18 leaves the Mondo term kg2 lists it on, viral hepatitis in animals, for its twin
    assert by_id["ICD10:B18"] is by_id["ICD10CM:B18"]
    assert by_id["MONDO:0025085"]["equivalent_ids"] == ["MONDO:0025085"]
    # R10's chapter says symptom, but it is its ICD-10-CM code: the node's type is that code's. A code that did not
    # move keeps its own
    assert by_id["ICD10:R10"]["categories"] == [DISEASE]
    assert by_id["ICD10:Y06.0"]["categories"] == ["biolink:Event"]
    # A code some source gives a species is not moved: it would give its twin's node that species
    assert by_id["ICD10:G61.0"] is not by_id["ICD10CM:G61.0"]
    # G50.1 joins Mondo's term through its ICD-10-CM twin. A name to match on would be evidence between the code and
    # the concept of that name, UMLS's, which Mondo's term is kept apart from
    facial_pain = by_id["MONDO:0018362"]
    assert by_id["ICD10:G50.1"] is facial_pain and by_id["ICD10CM:G50.1"] is facial_pain
    assert by_id["UMLS:C0154729"] is not facial_pain
    assert "Atypical facial pain" in facial_pain["synonyms"]
    # What Babel places stays as it was: E11 on type 2 diabetes, and Q93.5 on its syndrome without its ICD-10-CM
    # twin or its WHO name following it there
    assert {"ICD10:E11", "ICD10CM:E11"} <= set(by_id["MONDO:0005148"]["equivalent_ids"])
    syndrome = by_id["MONDO:0012315"]
    assert by_id["ICD10:Q93.5"] is syndrome and by_id["ICD10CM:Q93.5"] is not syndrome
    assert "Other deletions of part of a chromosome" not in syndrome.get("synonyms", [])
    # I21, named like the disease, takes the one place an ICD-10-CM code has in it. I21.9, which Mondo calls exact to
    # it and so has no name to match on, stays out with its WHO twin, and the two are named all the same
    infarction = by_id["MONDO:0004781"]
    assert by_id["ICD10CM:I21"] is infarction and by_id["ICD10:I21"] is infarction
    left_out = by_id["ICD10CM:I21.9"]
    assert left_out is not infarction and by_id["ICD10:I21.9"] is left_out
    assert left_out["name"] == "Acute myocardial infarction, unspecified"
    assert "Acute myocardial infarction NOS" in left_out["synonyms"], "kg2's name for the WHO code is kept"


def test_a_twin_statement_moves_only_the_who_code(tmp_path):
    """Every other id sits with the same ids as in a build without the statements, and the node a code joins keeps
    its id, name, types, description and taxon: here every node does. Weighed as evidence, a statement took its
    ICD-10-CM code out of the node its name had put it in."""
    with_statements = _resolve(tmp_path / "with")
    without = _resolve(tmp_path / "without", twin_statements=False)
    statements = jsonlines.open(tmp_path / "with" / "icd10" / "edges.jsonl")
    who_twins = {edge["subject"] for edge in statements if edge["predicate"] == SAME_AS}

    def nodes(by_id: dict[str, dict]) -> dict[frozenset[str], tuple]:
        return {
            frozenset(node["equivalent_ids"])
            - who_twins: (
                *(node["id"], node.get("name"), node["categories"]),
                *(node.get("description"), node.get("taxon")),
            )
            for node in by_id.values()
            if not who_twins.issuperset(node["equivalent_ids"])
        }

    assert nodes(with_statements) == nodes(without)
