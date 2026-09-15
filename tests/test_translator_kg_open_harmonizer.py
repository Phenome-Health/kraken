"""translator-kg-open nodes: equivalent_identifiers become equivalent_ids, and the SMILES in `xref` join them."""

from kraken.harmonizers.translator_kg_open import TranslatorKGOpenHarmonizer
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
