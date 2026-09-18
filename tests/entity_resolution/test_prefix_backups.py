"""Last-resort category/taxon guesses from an identifier's prefix."""

from kraken.entity_resolution.prefix_backups import infer_category, infer_taxon


def test_infer_category():
    assert infer_category("HGNC:2707") == "biolink:Gene"
    assert infer_category("UniProtKB:P12821") == "biolink:Protein"
    assert infer_category("CHEBI:1234") == "biolink:ChemicalEntity"
    assert infer_category("REACT:R-HSA-1") is None  # spans pathways, reactions and physical entities
    assert infer_category("WEIRD:1") is None


def test_infer_taxon_only_for_single_species_authorities():
    assert infer_taxon("MGI:98834") == "NCBITaxon:10090"
    assert infer_taxon("HGNC:2707") == "NCBITaxon:9606"
    assert infer_taxon("NCBIGene:1636") is None  # spans every species
