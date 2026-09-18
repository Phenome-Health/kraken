"""Tests for TaxonNormalizer -- species rollup and scientific-name lookup over NCBI's taxdump."""

import pytest

from kraken.utils.taxonomy import TaxonNormalizer
from tests.helpers import build_test_taxonomy, write_test_taxdump

# A small slice of the real tree, covering each case the rollup has to handle:
#   559292 (strain)      -> 4932   -- the case that motivates all of this
#   511145 (strain)      -> 562    -- two levels below its species
#   9606   (species)     -> itself
#   2      (superkingdom)-> itself -- above species; nothing to roll up to
#   10239  (no rank)     -> itself -- viruses; no species ancestor
NODES = [
    ("1", "1", "no rank"),
    ("2", "1", "superkingdom"),
    ("4932", "1", "species"),
    ("559292", "4932", "strain"),
    ("562", "2", "species"),
    ("83333", "562", "strain"),
    ("511145", "83333", "strain"),
    ("9606", "1", "species"),
    ("10239", "1", "no rank"),
]
NAMES = {
    "4932": "Saccharomyces cerevisiae",
    "559292": "Saccharomyces cerevisiae S288C",
    "562": "Escherichia coli",
    "9606": "Homo sapiens",
}


@pytest.fixture
def taxonomy(tmp_path) -> TaxonNormalizer:
    return build_test_taxonomy(tmp_path, NODES, NAMES)


def test_strain_rolls_up_to_species(taxonomy):
    assert taxonomy.to_species("559292") == "4932"


def test_rollup_climbs_multiple_levels(taxonomy):
    """E. coli K-12 MG1655 sits two strain levels below the species."""
    assert taxonomy.to_species("511145") == "562"
    assert taxonomy.to_species("83333") == "562"


def test_species_maps_to_itself(taxonomy):
    assert taxonomy.to_species("9606") == "9606"


def test_taxa_above_species_are_unchanged(taxonomy):
    """A genus/superkingdom can't be pushed down to a species, so it's left alone."""
    assert taxonomy.to_species("2") == "2"


def test_taxa_with_no_species_ancestor_are_unchanged(taxonomy):
    assert taxonomy.to_species("10239") == "10239"


def test_unknown_taxon_is_unchanged(taxonomy):
    assert taxonomy.to_species("999999") == "999999"


def test_curie_form_is_preserved(taxonomy):
    """Callers may pass a bare id or a CURIE; they get back the form they gave."""
    assert taxonomy.to_species("NCBITaxon:559292") == "NCBITaxon:4932"
    assert taxonomy.to_species("NCBITaxon:9606") == "NCBITaxon:9606"


def test_scientific_name_lookup(taxonomy):
    assert taxonomy.scientific_name("559292") == "Saccharomyces cerevisiae S288C"
    assert taxonomy.scientific_name("NCBITaxon:9606") == "Homo sapiens"


def test_scientific_name_skips_other_name_classes(taxonomy):
    """names.dmp carries synonyms and common names too; only the scientific name is kept."""
    assert taxonomy.scientific_name("9606") == "Homo sapiens"


def test_scientific_name_unknown_returns_none(taxonomy):
    assert taxonomy.scientific_name("999999") is None


def test_only_sub_species_taxa_are_retained(taxonomy):
    """Species and above map to themselves, so keeping them would waste memory across ~2.7M NCBI taxa."""
    assert taxonomy._species_of == {"559292": "4932", "83333": "562", "511145": "562"}


def test_missing_taxdump_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        TaxonNormalizer(tmp_path / "nope.tar.gz")


def test_cyclic_dump_does_not_hang(tmp_path):
    """A malformed dump shouldn't send the rollup into an infinite climb."""
    taxdump = write_test_taxdump(tmp_path, [("100", "200", "strain"), ("200", "100", "strain")], {})
    assert TaxonNormalizer(taxdump).to_species("100") == "100"
