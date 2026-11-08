"""translator-kg-open nodes: equivalent_identifiers become equivalent_ids, and the SMILES in `xref` join them."""

import json

import pytest

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.translator_kg_open import TranslatorKGOpenHarmonizer, release_babel_version
from tests.helpers import stub_normalization


class _LeafBiolink:
    def filter_to_leaf_categories(self, categories):
        return list(categories)


def _harmonizer() -> TranslatorKGOpenHarmonizer:
    """A harmonizer allocated without __init__ (which builds a Biolink toolkit and biomapper2 Normalizer)."""
    harmonizer = object.__new__(TranslatorKGOpenHarmonizer)
    harmonizer.source_infores = "translator-kg-open"
    harmonizer.biolink = _LeafBiolink()
    harmonizer.name_override_count = 0
    harmonizer.multi_taxon_node_count = 0
    harmonizer.multi_taxon_examples = []
    harmonizer.core_node_props = (
        {harmonizer.id_prop, harmonizer.category_prop, harmonizer.equivalent_ids_prop, harmonizer.name_prop}
        | harmonizer.synonyms_props
        | harmonizer.taxon_props
    )
    return stub_normalization(harmonizer)


def test_equivalent_identifiers_and_xref_smiles_become_equivalent_ids():
    node = _harmonizer()._harmonize_node(
        {
            "id": "CHEBI:100147",
            "category": ["biolink:SmallMolecule"],
            "name": "nalidixic acid",
            "equivalent_identifiers": ["CHEBI:100147", "PUBCHEM.COMPOUND:4421"],
            "xref": [
                "SMILES:CCn1cc(C(=O)O)c(=O)c2ccc(C)nc21",
                "InChIKey:MHWLWQUZZRMNGJ-UHFFFAOYSA-N",  # a case variant of an id already listed: not folded in
                "CHEMBL.TARGET:CHEMBL1234",  # a target, not this compound: not folded in
            ],
        }
    )
    assert set(node["equivalent_ids"]) == {
        "CHEBI:100147",
        "PUBCHEM.COMPOUND:4421",
        "SMILES:CCn1cc(C(=O)O)c(=O)c2ccc(C)nc21",
    }
    assert "xref" in node["attributes"]["translator-kg-open"]  # the whole xref is still kept


def test_node_without_smiles_keeps_just_its_equivalent_identifiers():
    node = _harmonizer()._harmonize_node(
        {"id": "MONDO:1", "category": ["biolink:Disease"], "equivalent_identifiers": ["MONDO:1", "DOID:2"]}
    )
    assert set(node["equivalent_ids"]) == {"MONDO:1", "DOID:2"}


def test_babel_version_is_read_from_the_metadata_beside_the_nodes(tmp_path):
    (tmp_path / "release-metadata.json").write_text(json.dumps({"babel_version": "2026jul22"}))
    assert release_babel_version(tmp_path / "nodes.jsonl") == "2026jul22"


@pytest.mark.parametrize(
    "metadata",
    [{}, {"babel_version": None}, {"babel_version": ""}, {"babel_version": "2026jul22 "}, {"babel_version": 20260722}],
)
def test_a_release_that_does_not_name_its_babel_is_an_error_not_a_guess(tmp_path, metadata):
    (tmp_path / "release-metadata.json").write_text(json.dumps(metadata))
    with pytest.raises(ValueError, match="babel_version"):
        release_babel_version(tmp_path / "nodes.jsonl")


def test_a_release_without_metadata_is_an_error_not_a_guess(tmp_path):
    with pytest.raises(FileNotFoundError):
        release_babel_version(tmp_path / "nodes.jsonl")


def test_harmonize_puts_the_release_metadata_beside_the_harmonized_files(tmp_path, monkeypatch):
    """Entity resolution reads only harmonized files, so the Babel release goes with them. A release that doesn't name
    its Babel fails before the harmonize, not after it."""
    harmonized = []
    monkeypatch.setattr(BaseHarmonizer, "harmonize", lambda self, nodes_output, *_a, **_k: harmonized.append(1))
    raw, out = tmp_path / "raw", tmp_path / "harmonized"
    raw.mkdir()
    out.mkdir()
    inputs = {"nodes_input": raw / "nodes.jsonl", "edges_input": raw / "edges.jsonl"}

    (raw / "release-metadata.json").write_text(json.dumps({"babel_version": "2026jul22"}))
    _harmonizer().harmonize(out / "nodes.jsonl", out / "edges.jsonl", **inputs)
    assert harmonized and release_babel_version(out / "nodes.jsonl") == "2026jul22"

    harmonized.clear()
    (raw / "release-metadata.json").write_text(json.dumps({}))
    with pytest.raises(ValueError, match="babel_version"):
        _harmonizer().harmonize(out / "nodes.jsonl", out / "edges.jsonl", **inputs)
    assert not harmonized
