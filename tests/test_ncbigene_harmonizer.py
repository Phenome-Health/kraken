"""Tests for the NCBI Gene harmonizer's file discovery and TSV parsing.

These deliberately avoid constructing a real harmonizer: BaseHarmonizer.__init__ builds a Biolink toolkit and
a biomapper2 Normalizer, both of which reach the network. Everything exercised here is independent of that, so
the instance is allocated without running __init__ (see `harmonizer`). The row -> node conversion, which does
need Biolink, is covered by running a real build against NCBI's files rather than in unit tests.
"""

import gzip
from collections import defaultdict

import pytest

from kraken.harmonizers.ncbigene import NCBIGeneHarmonizer
from kraken.harmonizers.ncbigene_taxon_allowlist import TAXON_ALLOWLIST
from tests.helpers import build_test_taxonomy

GENE_INFO_HEADER = (
    "#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\t"
    "type_of_gene\tSymbol_from_nomenclature_authority\tFull_name_from_nomenclature_authority\t"
    "Nomenclature_status\tOther_designations\tModification_date\tFeature_type"
)
APP_ROW = (
    "9606\t351\tAPP\t-\tAAA|ABETA\tMIM:104760|HGNC:HGNC:620|Ensembl:ENSG00000142192|AllianceGenome:HGNC:620\t"
    "21\t21q21.3\tamyloid beta precursor protein\tprotein-coding\tAPP\tamyloid beta precursor protein\tO\t"
    "amyloid-beta precursor protein\t20260809\t-"
)


class _StubNormalizer:
    """Stands in for biomapper2's normalizer, whose real one needs a Biolink toolkit (network). Knows two
    prefixes and accepts only numeric ids -- enough to exercise how the harmonizer feeds it and reads its
    three return values. biomapper2's own vocab coverage is tested in that repo, not here."""

    CANONICAL_PREFIXES = {"HGNC": "HGNC", "MIM": "OMIM"}

    def get_curies(self, local_ids_dict, **kwargs):
        curies, invalid_ids, unrecognized = {}, {}, set()
        for prefix, local_ids in local_ids_dict.items():
            if prefix not in self.CANONICAL_PREFIXES:
                unrecognized.add(prefix)
                continue
            for local_id in local_ids:
                if local_id.isdigit():
                    curies[f"{self.CANONICAL_PREFIXES[prefix]}:{local_id}"] = ""
                else:
                    invalid_ids.setdefault(prefix, []).append(local_id)
        return curies, invalid_ids, unrecognized


@pytest.fixture
def harmonizer() -> NCBIGeneHarmonizer:
    """A harmonizer allocated without __init__ (which would hit the network), with just the counters the
    methods under test touch."""
    instance = object.__new__(NCBIGeneHarmonizer)
    instance.normalizer = _StubNormalizer()
    instance.unrecognized_dbxref_prefixes = defaultdict(int)
    instance.invalid_dbxref_ids = defaultdict(int)
    instance.skipped_placeholders = 0
    instance.skipped_biological_regions = 0
    instance.skipped_uncharacterized = 0
    instance.unrecognized_gene_types = defaultdict(int)
    instance.skipped_by_taxon = 0
    instance.taxonomy = None
    instance.rolled_up_to_species = 0
    instance.use_taxon_allowlist = False  # tests that care about it opt in explicitly
    return instance


def _taxonomy(directory):
    """A TaxonNormalizer over a tiny taxdump: S288C (strain) under S. cerevisiae (species), plus human/cow."""
    return build_test_taxonomy(
        directory,
        nodes=[("1", "1", "no rank"), ("4932", "1", "species"), ("559292", "4932", "strain"), ("9606", "1", "species")],
        names={"4932": "Saccharomyces cerevisiae", "559292": "Saccharomyces cerevisiae S288C", "9913": "Bos taurus"},
    )


def _write_gene_info(path, rows=(APP_ROW,), gzipped=False):
    content = "\n".join([GENE_INFO_HEADER, *rows]) + "\n"
    if gzipped:
        path.write_bytes(gzip.compress(content.encode()))
    else:
        path.write_text(content)
    return path


# ------------------------------- input discovery -------------------------------


def test_finds_all_species_file(harmonizer, tmp_path):
    _write_gene_info(tmp_path / "gene_info")
    assert harmonizer._find_gene_info_files(tmp_path) == [tmp_path / "gene_info"]


def test_finds_gzipped_and_per_group_files(harmonizer, tmp_path):
    _write_gene_info(tmp_path / "All_Mammalia.gene_info.gz", gzipped=True)
    _write_gene_info(tmp_path / "Fungi.gene_info.gz", gzipped=True)
    found = harmonizer._find_gene_info_files(tmp_path)
    assert [p.name for p in found] == ["All_Mammalia.gene_info.gz", "Fungi.gene_info.gz"]


def test_rejects_mixing_all_species_with_per_group_files(harmonizer, tmp_path):
    """The per-group files are subsets of the all-species file, so ingesting both would double up genes."""
    _write_gene_info(tmp_path / "gene_info.gz", gzipped=True)
    _write_gene_info(tmp_path / "All_Mammalia.gene_info.gz", gzipped=True)
    with pytest.raises(ValueError, match="both the all-species"):
        harmonizer._find_gene_info_files(tmp_path)


def test_raises_when_no_gene_info_present(harmonizer, tmp_path):
    (tmp_path / "gene_summary.gz").write_bytes(gzip.compress(b"#tax_id\tGeneID\tSource\tSummary\n"))
    with pytest.raises(FileNotFoundError, match="no gene_info file found"):
        harmonizer._find_gene_info_files(tmp_path)


def test_prefers_uncompressed_over_gzipped(harmonizer, tmp_path):
    """An already-unzipped copy avoids decompressing several GB on every build."""
    _write_gene_info(tmp_path / "gene_info")
    _write_gene_info(tmp_path / "gene_info.gz", gzipped=True)
    assert harmonizer._resolve_optional(tmp_path, "gene_info") == tmp_path / "gene_info"


# ------------------------------- TSV parsing -------------------------------


def test_streams_tsv_stripping_commented_header(harmonizer, tmp_path):
    """NCBI comments out the header line, so '#tax_id' must become a usable 'tax_id' column."""
    rows = list(harmonizer._stream_tsv(_write_gene_info(tmp_path / "gene_info")))
    assert len(rows) == 1
    assert rows[0]["tax_id"] == "9606"
    assert rows[0]["Symbol"] == "APP"


def test_streams_gzipped_tsv_without_unzipping(harmonizer, tmp_path):
    rows = list(harmonizer._stream_tsv(_write_gene_info(tmp_path / "gene_info.gz", gzipped=True)))
    assert [row["GeneID"] for row in rows] == ["351"]


def test_skips_malformed_rows_rather_than_misaligning_columns(harmonizer, tmp_path):
    path = _write_gene_info(tmp_path / "gene_info", rows=(APP_ROW, "9606\t352\ttruncated"))
    assert [row["GeneID"] for row in harmonizer._stream_tsv(path)] == ["351"]


def test_value_normalizes_ncbi_null_marker(harmonizer):
    row = {"LocusTag": "-", "Symbol": "APP"}
    assert harmonizer._value(row, "LocusTag") == ""
    assert harmonizer._value(row, "Symbol") == "APP"
    assert harmonizer._value(row, "missing_column") == ""


def test_values_splits_pipe_delimited_lists(harmonizer):
    assert harmonizer._values({"Synonyms": "AAA|ABETA|AD1"}, "Synonyms") == ["AAA", "ABETA", "AD1"]
    assert harmonizer._values({"Synonyms": "-"}, "Synonyms") == []


# ------------------------------- dbXrefs -------------------------------


def test_equivalent_ids_delegates_to_biomapper2_and_fixes_doubled_prefixes(harmonizer):
    """Prefix/id validity is biomapper2's call, so this checks the integration: NCBI's doubled prefixes are
    stripped before handing ids over, and whatever comes back is used verbatim (note MIM -> OMIM is
    biomapper2 canonicalizing, not this harmonizer remapping)."""
    row = {"dbXrefs": "MIM:104760|HGNC:HGNC:620"}
    assert sorted(harmonizer._equivalent_ids(row)) == ["HGNC:620", "OMIM:104760"]


def test_equivalent_ids_drops_ignored_prefixes_without_flagging_them(harmonizer):
    """AllianceGenome re-points at ids we already take, so it never reaches biomapper2 and isn't reported as
    an unrecognized vocab -- it's a deliberate exclusion, not a gap."""
    assert harmonizer._equivalent_ids({"dbXrefs": "AllianceGenome:HGNC:620"}) == []
    assert harmonizer.unrecognized_dbxref_prefixes == {}


def test_equivalent_ids_counts_prefixes_biomapper2_does_not_know(harmonizer):
    """Counted so a missing vocab surfaces in the run log as a one-line follow-up in biomapper2, rather
    than the xrefs vanishing."""
    assert harmonizer._equivalent_ids({"dbXrefs": "SomeNewDB:123|SomeNewDB:456"}) == []
    assert harmonizer.unrecognized_dbxref_prefixes == {"SomeNewDB": 1}


def test_equivalent_ids_counts_ids_that_fail_validation(harmonizer):
    """A known prefix carrying a malformed id is a data problem worth reporting, not a silent drop."""
    assert harmonizer._equivalent_ids({"dbXrefs": "HGNC:not-a-number"}) == []
    assert harmonizer.invalid_dbxref_ids == {"HGNC": 1}


# ------------------------------- uncharacterized genes -------------------------------
#
# Dropped only when a gene fails ALL THREE tests at once: no real symbol, no dbXref, no real description.


def test_uncharacterized_gene_is_dropped(harmonizer):
    assert harmonizer._is_uncharacterized("LOC38088273", "", "hypothetical protein", []) is True
    assert harmonizer._is_uncharacterized("B1U21_RS00005", "B1U21_RS00005", "", []) is True
    assert harmonizer._is_uncharacterized("LOC123", "", "uncharacterized protein", []) is True


def test_an_xref_alone_keeps_a_gene(harmonizer):
    """Linkable, so it can still resolve against other sources however poorly it's named."""
    assert harmonizer._is_uncharacterized("LOC38088273", "", "hypothetical protein", ["HGNC:620"]) is False


def test_a_real_symbol_alone_keeps_a_gene(harmonizer):
    """Findable by name -- and a real symbol is what other sources will match on."""
    assert harmonizer._is_uncharacterized("APP", "B1U21_RS00005", "hypothetical protein", []) is False


def test_a_real_description_alone_keeps_a_gene(harmonizer):
    """The case that rules out the blunter filter: rRNA/tRNA genes are unnamed and unlinked but real."""
    assert harmonizer._is_uncharacterized("LOC38088273", "", "5S ribosomal RNA", []) is False
    assert harmonizer._is_uncharacterized("OrniCt037", "OrniCt037", "tRNA-Ile", []) is False


def test_generic_description_matching_is_case_insensitive_and_prefix_based(harmonizer):
    """NCBI's phrasings vary ("Uncharacterized protein LOC123"), so matching is a lowercase prefix test."""
    assert harmonizer._is_uncharacterized("LOC123", "", "Uncharacterized protein LOC123", []) is True
    assert harmonizer._is_uncharacterized("LOC123", "", "Hypothetical protein, unlikely", []) is True


# ------------------------------- categories / filtering -------------------------------


@pytest.mark.parametrize("gene_type", ["protein-coding", "ncRNA", "pseudo", "tRNA", "unknown", "other"])
def test_gene_types_map_to_gene(harmonizer, gene_type):
    assert harmonizer._categories_for(gene_type, []) == ["biolink:Gene"]


def test_unrecognized_gene_type_still_maps_to_gene_but_is_counted(harmonizer):
    """New type_of_gene values appear in NCBI releases; dropping those genes would be worse than typing them
    as Gene -- but it's counted so the guess surfaces in the run log rather than passing unnoticed."""
    assert harmonizer._categories_for("some-future-type", []) == ["biolink:Gene"]
    assert harmonizer.unrecognized_gene_types == {"some-future-type": 1}


def test_pseudogenes_are_typed_as_gene(harmonizer):
    """Biolink 4.2.5 has no Pseudogene class (and Gene has no descendants), and ROBOKOP types real
    pseudogenes as biolink:Gene -- so anything else would fragment merging."""
    assert harmonizer._categories_for("pseudo", []) == ["biolink:Gene"]
    assert harmonizer.unrecognized_gene_types == {}


def test_biological_regions_are_filtered_out_by_default(harmonizer):
    assert harmonizer._categories_for("biological-region", ["regulatory:enhancer"]) is None


def test_biological_regions_categorized_when_included(harmonizer):
    harmonizer.include_biological_regions = True
    assert harmonizer._categories_for("biological-region", ["regulatory:enhancer"]) == ["biolink:RegulatoryRegion"]
    assert harmonizer._categories_for("biological-region", ["misc_recomb:non_allelic_homologous"]) == [
        "biolink:BiologicalEntity"
    ]


# ------------------------------- summaries -------------------------------


def test_load_summaries_interns_duplicate_texts(harmonizer, tmp_path):
    """3.1M summary rows share only ~150k distinct texts, so identical strings must be shared, not copied."""
    boilerplate = "microRNAs (miRNAs) are short non-coding RNAs."
    (tmp_path / "gene_summary").write_text(
        "#tax_id\tGeneID\tSource\tSummary\n"
        f"9606\t1\tRefSeq\t{boilerplate}\n"
        f"9606\t2\tRefSeq\t{boilerplate}\n"
        "9606\t3\tRefSeq\tA distinct summary.\n"
    )
    summaries = harmonizer._load_summaries(tmp_path)
    assert summaries["1"] == boilerplate
    assert summaries["3"] == "A distinct summary."
    assert summaries["1"] is summaries["2"]  # interned, not two copies of the same text


def test_load_summaries_returns_empty_when_file_absent(harmonizer, tmp_path):
    harmonizer.use_full_name_as_description_fallback = True
    assert harmonizer._load_summaries(tmp_path) == {}


# ------------------------------- taxon allowlist -------------------------------


def test_allowlist_covers_the_core_model_organisms():
    """The curated list exists because no ranking over gene_info surfaces these -- ranking by gene count puts
    a corroboree frog second and omits rat entirely."""
    assert {"9606", "10090", "10116", "7955", "7227", "6239", "4932", "3702"}.issubset(TAXON_ALLOWLIST)


def test_allowlist_entries_are_all_plain_tax_ids():
    """Bare NCBI tax_ids, no 'NCBITaxon:' prefix -- they're compared against gene_info's tax_id column."""
    assert all(entry.isdigit() for entry in TAXON_ALLOWLIST)


def test_allowlist_may_hold_sub_species_entries(harmonizer, tmp_path):
    """Entries are written at whatever rank names the organism -- dog is published as the subspecies Canis
    lupus familiaris -- and the harmonizer rolls them up before comparing, so no entry must be species-rank."""
    assert "9615" in TAXON_ALLOWLIST  # Canis lupus familiaris, a subspecies of Canis lupus


def test_taxonomy_absent_when_no_taxdump_and_allowlist_off(harmonizer, tmp_path):
    assert harmonizer._load_taxonomy(tmp_path) is None


def test_missing_taxdump_raises_when_allowlist_is_on(harmonizer, tmp_path):
    """A species-level allowlist without the taxonomy would silently drop every strain-published gene."""
    harmonizer.use_taxon_allowlist = True
    with pytest.raises(FileNotFoundError, match="required"):
        harmonizer._load_taxonomy(tmp_path)


def test_strain_published_gene_kept_when_its_species_is_allowlisted(harmonizer, tmp_path):
    """Yeast's genes are published under 559292 (S288C) while 4932 is what's curated, so a raw comparison
    would drop essentially all of them."""
    harmonizer.taxonomy = _taxonomy(tmp_path)
    harmonizer.allowed_species = frozenset({"4932"})

    # The gate is "roll up, then test membership" -- asserted directly, since going further through
    # _harmonize_row would need the Biolink client this offline fixture deliberately omits.
    assert harmonizer._species_tax_id("559292") in harmonizer.allowed_species


def test_taxa_rolled_up_to_species(harmonizer, tmp_path):
    """NCBI Gene reports many microbes at strain rank; other sources report the species. Storing the strain
    unchanged would make a species-level `taxa` query silently miss these genes."""
    harmonizer.taxonomy = _taxonomy(tmp_path)
    assert harmonizer._species_tax_id("559292") == "4932"  # S288C -> S. cerevisiae
    assert harmonizer.rolled_up_to_species == 1


def test_species_and_unknown_taxa_pass_through(harmonizer, tmp_path):
    harmonizer.taxonomy = _taxonomy(tmp_path)
    assert harmonizer._species_tax_id("9606") == "9606"  # already a species
    assert harmonizer._species_tax_id("999999") == "999999"  # unknown to the taxonomy
    assert harmonizer._species_tax_id("") == ""
    assert harmonizer.rolled_up_to_species == 0


def test_taxa_unchanged_without_taxonomy(harmonizer):
    """No taxdump means no rollup -- the strain id is used as NCBI reported it."""
    assert harmonizer._species_tax_id("559292") == "559292"
    assert harmonizer.rolled_up_to_species == 0


def test_fallback_description_names_the_organism(harmonizer, tmp_path):
    """The same gene name recurs across species, so an unqualified full name reads as the human one."""
    harmonizer.taxonomy = _taxonomy(tmp_path)
    assert (
        harmonizer._full_name_description("amyloid beta precursor protein", "9913")
        == "amyloid beta precursor protein (Bos taurus)"
    )


def test_fallback_description_uses_the_reported_taxon_not_the_rolled_up_one(harmonizer, tmp_path):
    """A strain-specific gene should still say so in prose, even though `taxa` is rolled up to species."""
    harmonizer.taxonomy = _taxonomy(tmp_path)
    assert harmonizer._full_name_description("some gene", "559292") == "some gene (Saccharomyces cerevisiae S288C)"


def test_full_name_drops_descriptions_that_just_restate_the_gene_type(harmonizer):
    """~0.9% of genes have a `description` of just "pseudo" or "ncRNA" -- a category, not a name."""
    assert harmonizer._gene_full_name({"description": "pseudo"}, "pseudo") == ""
    assert harmonizer._gene_full_name({"description": "ncRNA"}, "ncRNA") == ""
    assert (
        harmonizer._gene_full_name({"description": "amyloid beta precursor protein"}, "protein-coding")
        == "amyloid beta precursor protein"
    )


def test_fallback_description_omits_parenthetical_for_unknown_taxon(harmonizer, tmp_path):
    harmonizer.taxonomy = _taxonomy(tmp_path)
    assert harmonizer._full_name_description("some gene", "999999") == "some gene"
    assert harmonizer._full_name_description("", "9913") == ""


def test_rows_outside_allowlist_are_filtered_out(harmonizer, tmp_path):
    harmonizer.taxonomy = _taxonomy(tmp_path)
    harmonizer.use_taxon_allowlist = True
    harmonizer.allowed_species = frozenset({"4932"})
    row = dict(zip(GENE_INFO_HEADER.lstrip("#").split("\t"), APP_ROW.split("\t")))

    assert harmonizer._harmonize_row({**row, "tax_id": "9606"}, {}) is None
    assert harmonizer.skipped_by_taxon == 1
