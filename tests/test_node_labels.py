"""Post-build node labels (post_processing/node_labels.py): what gets each label, and that nothing else changes."""

import json

import pytest

from kraken.post_processing.node_labels import label_nodes, template_family


def _pathway(name: str, curie: str = "PathWhiz:PW1") -> dict:
    return {"id": curie, "name": name, "categories": ["biolink:Pathway"], "equivalent_ids": [curie]}


CL, PC, TG = "Cardiolipin Biosynthesis", "Phosphatidylcholine Biosynthesis", "De Novo Triacylglycerol Biosynthesis"


@pytest.mark.parametrize(
    "name, family",
    [
        (f"{CL} CL(10:0/10:0/10:0/28:0) pathway", (CL, "CL(10:0/10:0/10:0/28:0)")),
        (f"{CL} CL(i-19:0/10:0(3-OH)/10:0/10:0) pathway", (CL, "CL(i-19:0/10:0(3-OH)/10:0/10:0)")),
        (f"{PC} PC(14:0/14:0) | PE(14:0/14:0)", (PC, "PC(14:0/14:0) | PE(14:0/14:0)")),
        (f"{TG} TG(10:0/10:0/10:0),Metabolic pathway", (TG, "TG(10:0/10:0/10:0),Metabolic")),
        ("Glycolysis pathway", None),  # no lipid in the name
        ("PC(14:0/14:0) pathway", None),  # nothing before the lipid to name a family
    ],
)
def test_template_family_from_the_name(name, family):
    assert template_family(_pathway(name)) == family


def test_only_pathwhiz_and_smpdb_pathways_are_templates():
    assert template_family(_pathway("Cardiolipin Biosynthesis CL(10:0/10:0)", "REACT:R-1")) is None
    # a PathWhiz compound id is not a PathWhiz pathway id, even on a node typed Pathway
    assert template_family(_pathway("Cardiolipin Biosynthesis CL(10:0/10:0)", "PathWhiz.Compound:1")) is None


def _write(path, rows):
    path.write_text("".join(json.dumps(row, ensure_ascii=False) + "\n" for row in rows))


def _read(path):
    return {row["id"]: row for row in map(json.loads, path.read_text().splitlines())}


def _chemical(curie, name):
    return {"id": curie, "name": name, "categories": ["biolink:SmallMolecule"], "equivalent_ids": [curie]}


def _vocabulary(curie, name):
    return {"id": curie, "name": name, "categories": ["biolink:ClinicalMeasurement"], "equivalent_ids": [curie]}


def _edge(subject, predicate, object_):
    return {"subject": subject, "predicate": f"biolink:{predicate}", "object": object_}


@pytest.fixture
def graph(tmp_path):
    nodes, edges = tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl"
    _write(
        nodes,
        [
            {**_chemical("HMDB:HMDB0035409", "Ozone"), "attributes": {"infores:x": {"kept": True}}},
            # a lab test with only a mapping edge: labelled, and its other attributes stay
            {**_vocabulary("LOINC:105571-4", "Ozone"), "attributes": {"infores:loinc": {"kept": True}}},
            _chemical("HMDB:HMDB0000062", "Carnitine"),
            _vocabulary("LOINC:LP14482-1", "Carnitine"),  # has biology of its own: not labelled
            _chemical("CHEBI:1", "Lonely"),
            _vocabulary("UMLS:C1", "Lonely"),  # its chemical has no biology: not labelled
            _pathway("Cardiolipin Biosynthesis CL(10:0/10:0/10:0/28:0) pathway", "PathWhiz:PW003484"),
            _pathway("Cardiolipin Biosynthesis CL(10:0/10:0/10:0/30:0) pathway", "PathWhiz:PW003501"),
            {"id": "MONDO:1", "name": "Ozone", "categories": ["biolink:Disease"], "equivalent_ids": ["MONDO:1"]},
            # a chemical node that also holds a vocabulary id is not vocabulary-only: not labelled, but a candidate
            {**_chemical("CHEBI:2", "Ozone"), "equivalent_ids": ["CHEBI:2", "UMLS:C2"]},
            # a vocabulary node with no edges at all is labelled, with every candidate in order
            _vocabulary("NCIT:C1", "Glycine"),
            _chemical("RM:2", "Glycine"),
            {**_chemical("UMLS:C3", "Glycine"), "equivalent_ids": ["UMLS:C3", "HMDB:HMDB0000123"]},
            # a chemical with a vocabulary id and no edges either: not vocabulary-only, so still not labelled
            {**_chemical("CHEBI:3", "Glycine"), "equivalent_ids": ["CHEBI:3", "UMLS:C4"]},
            _chemical("CHEBI:16236", "Ethanol"),
            _vocabulary("UMLS:C5", "Ethanol"),  # labelled: the control for the three below
            _vocabulary("LOINC:LP3", "Ethanol"),  # a biological edge pointing AT it: not labelled
            _vocabulary("NCIT:C2", "Ethanol"),  # an edge of a kind neither listed: not labelled
            {**_vocabulary("MESH:D1", "Ethanol"), "equivalent_ids": ["MESH:D1", "CHEBI:5"]},  # holds a chemical id
        ],
    )
    _write(
        edges,
        [
            _edge("HMDB:HMDB0035409", "participates_in", "PathWhiz:PW003484"),
            _edge("LOINC:105571-4", "related_to", "HMDB:HMDB0035409"),
            _edge("HMDB:HMDB0000062", "participates_in", "PathWhiz:PW003501"),
            _edge("LOINC:LP14482-1", "affects", "MONDO:1"),
            _edge("RM:2", "participates_in", "PathWhiz:PW003484"),
            _edge("UMLS:C3", "participates_in", "PathWhiz:PW003484"),
            _edge("CHEBI:2", "participates_in", "PathWhiz:PW003484"),
            _edge("CHEBI:16236", "participates_in", "PathWhiz:PW003484"),
            _edge("MONDO:1", "affects", "LOINC:LP3"),
            _edge("NCIT:C2", "treats", "MONDO:1"),
        ],
    )
    return nodes, edges


def test_labels(graph):
    nodes, edges = graph
    label_nodes(nodes, edges)
    result = _read(nodes)
    labels = {curie: node.get("attributes", {}).get("kraken") for curie, node in result.items()}
    assert labels["LOINC:105571-4"] == {"same_name_chemical_candidates": ["CHEBI:2", "HMDB:HMDB0035409"]}
    assert labels["NCIT:C1"] == {"same_name_chemical_candidates": ["RM:2", "UMLS:C3"]}
    assert labels["UMLS:C5"] == {"same_name_chemical_candidates": ["CHEBI:16236"]}
    assert result["LOINC:105571-4"]["attributes"]["infores:loinc"] == {"kept": True}
    for unlabelled in (
        "LOINC:LP14482-1",
        "UMLS:C1",
        "MONDO:1",
        "HMDB:HMDB0035409",
        "CHEBI:2",
        "UMLS:C3",
        "CHEBI:3",
        "LOINC:LP3",
        "NCIT:C2",
        "MESH:D1",
    ):
        assert "kraken" not in result[unlabelled].get("attributes", {})
    assert result["PathWhiz:PW003501"]["attributes"]["kraken"]["pathway_template_family"] == {
        "name": "Cardiolipin Biosynthesis",
        "variant": "CL(10:0/10:0/10:0/30:0)",
        "members": 2,
    }
    assert result["HMDB:HMDB0035409"]["attributes"] == {"infores:x": {"kept": True}}


def test_a_second_run_changes_nothing_and_a_stale_label_goes(graph):
    nodes, edges = graph
    label_nodes(nodes, edges)
    before = nodes.read_text()
    assert label_nodes(nodes, edges)["changed"] == 0
    assert nodes.read_text() == before
    _write(edges, [_edge("LOINC:105571-4", "affects", "MONDO:1")])  # now it has an edge of its own
    label_nodes(nodes, edges)
    assert _read(nodes)["LOINC:105571-4"]["attributes"] == {"infores:loinc": {"kept": True}}  # only the label went
