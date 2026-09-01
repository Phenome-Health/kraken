# ncbigene.py
import gzip
import logging
import re
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import jsonlines

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.ncbigene_taxon_allowlist import TAXON_ALLOWLIST
from kraken.utils.constants import NCBIGENE_INFORES
from kraken.utils.general import is_empty
from kraken.utils.kg_io import fix_repeated_prefix, save_to_jsonl
from kraken.utils.taxonomy import TaxonNormalizer

# --- Input file names (NCBI's own, as published under ftp.ncbi.nlm.nih.gov/gene/DATA/) ---
# The all-species bulk file. Mutually exclusive with the per-group files below (see _find_gene_info_files).
ALL_SPECIES_FILENAME = "gene_info"
# Per-taxonomic-group files from GENE_INFO/<group>/, e.g. "All_Mammalia.gene_info". Using these instead of
# the all-species file is how you scope the ingest (see the class docstring).
GROUP_FILE_SUFFIX = ".gene_info"
SUMMARY_FILENAME = "gene_summary"
# NCBI's taxonomy archive (ftp.ncbi.nlm.nih.gov/pub/taxonomy/), read straight out of the tarball: `nodes.dmp`
# to roll each gene's taxon up to species rank, `names.dmp` to name the organism in its description. The
# rollup matters because NCBI Gene publishes many microbes at STRAIN rank while other sources cite the species
# (NCBITaxon:559292 and :4932 are the same organism).
TAXDUMP_FILENAME = "taxdump.tar.gz"

# --- gene_info columns ---
COL_TAX_ID = "tax_id"  # published as "#tax_id"; the leading '#' is stripped when we read the header
COL_GENE_ID = "GeneID"
COL_SYMBOL = "Symbol"
COL_LOCUS_TAG = "LocusTag"
COL_SYNONYMS = "Synonyms"
COL_DBXREFS = "dbXrefs"
COL_CHROMOSOME = "chromosome"
COL_MAP_LOCATION = "map_location"
COL_DESCRIPTION = "description"
COL_TYPE_OF_GENE = "type_of_gene"
COL_AUTHORITY_SYMBOL = "Symbol_from_nomenclature_authority"
COL_AUTHORITY_FULL_NAME = "Full_name_from_nomenclature_authority"
COL_OTHER_DESIGNATIONS = "Other_designations"
COL_FEATURE_TYPE = "Feature_type"

# --- gene_summary columns ---
COL_SUMMARY = "Summary"

LIST_DELIMITER = "|"
EMPTY_VALUE = "-"  # NCBI's null marker in these files

# Placeholder records NCBI mints so GeneRIFs can be submitted for a gene that isn't in Gene yet. They carry no
# real gene identity (symbol is literally "NEWENTRY", no xrefs, no synonyms) -- roughly one per taxon.
PLACEHOLDER_SYMBOL = "NEWENTRY"

# --- type_of_gene -> Biolink category ---
# Every value here denotes an NCBI *gene* record, so all map to biolink:Gene -- including the RNA-gene types
# (a gene transcribed to ncRNA is still a Gene in Biolink; the RNA *product* would be a separate node).
#
# `pseudo` is the one uneasy fit: Biolink defines Gene as a region encoding a FUNCTIONAL transcript, which a
# pseudogene by definition does not. Biolink 4.2.5 offers nothing better though -- it has no Pseudogene class
# and Gene has no descendants -- and ROBOKOP types real pseudogenes (EEF1A1P49, SNRPGP16) as biolink:Gene, so
# doing otherwise would fragment merging for ~9% of these nodes. type_of_gene is kept in attributes, so this
# is revisitable if Biolink ever adds the class.
GENE_TYPES = {
    "protein-coding",
    "ncRNA",
    "snoRNA",
    "snRNA",
    "rRNA",
    "tRNA",
    "scRNA",
    "miscRNA",
    "pseudo",
    "unknown",
    "other",
}
GENE_CATEGORY = "biolink:Gene"

# `biological-region` records are NOT genes -- they're enhancers, silencers, imprinting control regions and
# similar annotated loci. They're excluded by default (see include_biological_regions): in the human file they
# are 66% of all rows yet only 16 of 128k carry any dbXref, so they add bulk and contribute nothing to entity
# resolution. Set include_biological_regions=True to ingest them under the categories below.
BIOLOGICAL_REGION_TYPE = "biological-region"
REGULATORY_REGION_CATEGORY = "biolink:RegulatoryRegion"
BIOLOGICAL_ENTITY_CATEGORY = "biolink:BiologicalEntity"
REGULATORY_FEATURE_PREFIX = "regulatory:"

# --- dbXrefs prefix handling ---
# Which prefixes are valid, and what their canonical form is, is biomapper2's job -- every xref goes through
# its normalizer so that knowledge lives in one place (biomapper2 core/normalizer/vocab_config.py) rather than
# being restated per harmonizer. Prefixes it doesn't know are counted and logged at the end of the run, so
# adding a vocab there is a visible, one-line follow-up rather than a silent loss.
#
# The only prefixes hardcoded here are ones that are semantically wrong as equivalences no matter what
# biomapper2 knows about them:
#   AllianceGenome -- not its own identifier space; it re-points at the HGNC/MGI/Xenbase ids we already take,
#                     so it would restate an equivalence we have (and its ids embed a second prefix)
#   AnimalQTLdb    -- quantitative trait loci, not gene identifiers; a QTL is not the same entity as a gene
DBXREF_PREFIXES_IGNORED = {
    "AllianceGenome",
    "AnimalQTLdb",
}

# --- Uncharacterized genes (see drop_uncharacterized_genes) ---
# When no nomenclature authority has named a gene, NCBI's Symbol is either the locus tag or a "LOC<GeneID>"
# placeholder. Neither is a real name, though both are at least unique -- which is why they're kept as the
# node's name rather than swapped for the description (see the class docstring).
LOC_PLACEHOLDER_PATTERN = re.compile(r"^LOC\d+$")
# Descriptions that state only that nothing is known. Matched as a lowercase prefix. Necessarily approximate,
# so the count of what it drops is logged each run and the list extended when new phrasings show up.
GENERIC_DESCRIPTION_PREFIXES = (
    "hypothetical protein",
    "uncharacterized",
    "putative uncharacterized",
    "unnamed protein",
    "unknown protein",
    "predicted protein",
    "protein of unknown function",
)

NCBI_GENE_PREFIX = "NCBIGene"
NCBI_TAXON_PREFIX = "NCBITaxon"
# Attribute holding the taxon as NCBI Gene reported it, recorded only when it differed from the species-rank
# value stored in `taxon` (i.e. for strain/subspecies records).
RAW_TAXON_ATTRIBUTE = "ncbi_reported_taxon"
GENE_URL_TEMPLATE = "https://www.ncbi.nlm.nih.gov/gene/{gene_id}"


class NCBIGeneHarmonizer(BaseHarmonizer):
    """Harmonizer for NCBI Gene's bulk files (ftp.ncbi.nlm.nih.gov/gene/DATA/).

    Ingests genes as nodes; emits no edges. Its value to KRAKEN is threefold:
      1. `tax_id` gives an authoritative taxon for every NCBIGene node. NCBIGene, unlike HGNC/MGI/RGD/ZFIN, is a
         species-agnostic prefix, so taxon cannot be inferred from the CURIE -- this file is the only way to get
         it in bulk, and it's what makes a one-taxon-per-entity rule bite on the largest gene population.
      2. `gene_summary` supplies real prose descriptions, which the aggregator KGs largely don't carry.
      3. `dbXrefs` is a clean, high-trust set of gene equivalencies (HGNC, Ensembl, OMIM, MGI, ...).

    Input is a DIRECTORY (config `input_file`) holding these files, gzipped or not:
      - `gene_info`      -- the all-species bulk file, OR one or more per-group `<group>.gene_info` files from
        GENE_INFO/<group>/ (e.g. `All_Mammalia.gene_info`). Mixing the all-species file with per-group files
        is rejected, since it would double-ingest genes.
      - `gene_summary`   -- optional; prose descriptions, joined on GeneID.
      - `taxdump.tar.gz` -- NCBI's taxonomy; required unless use_taxon_allowlist is off (see below).

    SCOPE. The all-species file is ~72M genes across ~54k taxa -- several times the size of the rest of KRAKEN
    combined, and overwhelmingly genomes of no biomedical interest. By default the ingest is restricted to the
    curated organisms in `ncbigene_taxon_allowlist.py` (see that module for why the list is curated rather
    than derived from any ranking). Set use_taxon_allowlist=False to take everything.

    Both the allowlist and each gene's taxon are rolled up to species before comparison, which is why the
    taxdump is required alongside it: NCBI publishes many genes under strain-level taxa (essentially all of
    yeast's live under 559292 / S288C, not 4932), so a raw comparison would silently drop them.

    Note that filtering on "has a dbXref" would NOT be a good proxy for relevance -- only 62.9% of the
    NCBIGene ids already in KRAKEN 2.1.0 carry one, so it would drop over a third of the very genes this
    harmonizer exists to annotate.

    NAMES. A gene's `name` is its symbol -- the nomenclature authority's if there is one, else NCBI's. For
    unnamed genes NCBI's symbol is a locus tag or a "LOC<GeneID>" placeholder, which looks like junk but is
    unique, and the functional text ("tRNA-Ile") is still carried in `description` and `synonyms`. Promoting
    that text to `name` was considered and rejected: it collides heavily (125k genes would be named
    "hypothetical protein"), and KG2 and ROBOKOP both name their NCBIGene nodes by symbol -- including ugly
    ones like EEF1A1P49 -- so a merged node's label would otherwise flip depending on merge order.

    Records filtered out regardless of scope: NEWENTRY placeholders, `biological-region` records (unless
    include_biological_regions=True), and fully uncharacterized genes (unless drop_uncharacterized_genes=False)
    -- see the constants and flags above for why.
    """

    source_infores = NCBIGENE_INFORES

    # Restrict the ingest to the curated organisms in ncbigene_taxon_allowlist.py (~2.4M genes; measure a
    # build's coverage with scripts/audit_ncbigene_taxon_allowlist.py). Set False to ingest EVERY taxon in the
    # input files, which for NCBI's all-species file means ~72M genes -- several times the rest of KRAKEN.
    use_taxon_allowlist: bool = True

    # Drop genes that are uncharacterized on every axis at once: no real symbol (NCBI fell back to a locus tag
    # or a LOC<GeneID> placeholder), no dbXref to link them by, and no description beyond a generic stand-in
    # like "hypothetical protein". Such a node can't be found by name, resolved by equivalence, or read -- it
    # is pure bulk. On the 2026-08-31 all-species file this dropped 535k of 2.28M genes (23%). Genes failing
    # only one or two of the three are KEPT: an rRNA gene with no xref still carries real content.
    drop_uncharacterized_genes: bool = True

    # Ingest `biological-region` records (enhancers, silencers, etc.). Off by default -- they are not genes,
    # and essentially none of them carry an xref, so they add bulk without helping entity resolution.
    include_biological_regions: bool = False

    # When a gene has no prose summary in gene_summary, fall back to gene_info's `description` column for the
    # node's description. NOTE that column is really a *full name* ("amyloid beta precursor protein"), not
    # prose -- so this trades field purity for coverage. Set False to leave such nodes without a description
    # (they still get the full name as a synonym either way).
    use_full_name_as_description_fallback: bool = True

    def __init__(self, biolink_client):
        super().__init__(biolink_client)
        self.unrecognized_dbxref_prefixes: dict[str, int] = defaultdict(int)
        self.invalid_dbxref_ids: dict[str, int] = defaultdict(int)
        self.skipped_placeholders = 0
        self.skipped_biological_regions = 0
        self.skipped_uncharacterized = 0
        self.unrecognized_gene_types: dict[str, int] = defaultdict(int)
        self.skipped_by_taxon = 0
        self.taxonomy: TaxonNormalizer | None = None  # None if no taxdump provided
        self.allowed_species: frozenset[str] = frozenset()  # TAXON_ALLOWLIST, rolled up to species rank
        self.rolled_up_to_species = 0

    def harmonize(
        self,
        nodes_output: Path,
        edges_output: Path,
        *,
        input_file: Path | None = None,
        nodes_input: Path | None = None,
        edges_input: Path | None = None,
    ):
        if not input_file:
            raise ValueError(f"{self.source_name} requires input_file (the directory holding NCBI's gene files)")
        input_dir = Path(input_file)
        if not input_dir.is_dir():
            raise ValueError(f"{self.source_name}: input_file must be a directory, got {input_dir}")

        gene_info_paths = self._find_gene_info_files(input_dir)
        logging.info(f"Harmonizing {self.source_name}: {[p.name for p in gene_info_paths]} -> {nodes_output}")

        self.taxonomy = self._load_taxonomy(input_dir)
        if self.use_taxon_allowlist:
            # Roll the allowlist up too, so its entries can be written at whatever rank names the organism
            # (dog is published as the subspecies Canis lupus familiaris) and still match a gene's species.
            self.allowed_species = frozenset(self.taxonomy.to_species(tax_id) for tax_id in TAXON_ALLOWLIST)
            logging.info(
                f"Ingesting the {len(TAXON_ALLOWLIST)} curated organisms in ncbigene_taxon_allowlist.py "
                f"({len(self.allowed_species)} distinct species after rollup)"
            )
        else:
            logging.info("Ingesting EVERY taxon in the input files (use_taxon_allowlist=False)")
        summaries = self._load_summaries(input_dir)

        node_count = 0
        with jsonlines.open(nodes_output, "w") as writer:
            for gene_info_path in gene_info_paths:
                for row in self._stream_tsv(gene_info_path):
                    node = self._harmonize_row(row, summaries)
                    if node:
                        writer.write(node)
                        node_count += 1

        save_to_jsonl([], edges_output, mode="w")  # NCBI Gene contributes nodes only

        self._log_run_summary(node_count)
        logging.info(f"{self.source_name} harmonization complete: {node_count} nodes, 0 edges")

    # ------------------------------- input discovery / streaming -------------------------------

    def _find_gene_info_files(self, input_dir: Path) -> list[Path]:
        """Locate the gene_info file(s) in the input dir, accepting gzipped or plain. Either the all-species
        `gene_info` file or one-or-more per-group `<group>.gene_info` files -- never both, since the groups are
        subsets of the all-species file and ingesting both would emit every grouped gene twice."""
        all_species = self._resolve_optional(input_dir, ALL_SPECIES_FILENAME)
        group_files = sorted(
            path
            for path in input_dir.iterdir()
            if path.name.removesuffix(".gz").endswith(GROUP_FILE_SUFFIX) and path != all_species
        )

        if all_species and group_files:
            raise ValueError(
                f"{self.source_name}: {input_dir} holds both the all-species '{ALL_SPECIES_FILENAME}' file and "
                f"per-group files ({[p.name for p in group_files]}). The groups are subsets of the all-species "
                f"file, so ingesting both would double up every grouped gene. Keep one or the other."
            )
        if all_species:
            return [all_species]
        if group_files:
            return group_files
        raise FileNotFoundError(
            f"{self.source_name}: no gene_info file found in {input_dir}. Expected '{ALL_SPECIES_FILENAME}' "
            f"(optionally .gz) or one or more per-group '*{GROUP_FILE_SUFFIX}' files."
        )

    @staticmethod
    def _resolve_optional(input_dir: Path, stem: str) -> Path | None:
        """Return the path to `stem` in input_dir, preferring an uncompressed copy over a .gz one (both are
        readable, but an already-unzipped file avoids decompressing on every build). None if neither exists."""
        for candidate in (input_dir / stem, input_dir / f"{stem}.gz"):
            if candidate.is_file():
                return candidate
        return None

    @staticmethod
    def _open_maybe_gzipped(path: Path):
        """Open a possibly-gzipped text file. We read NCBI's .gz files directly rather than unzipping them to
        disk first -- gene_info expands from ~1.5GB to several GB, and we only ever need one pass."""
        if path.suffix == ".gz":
            return gzip.open(path, "rt", encoding="utf-8", newline="")
        return open(path, encoding="utf-8", newline="")

    def _stream_tsv(self, path: Path) -> Iterator[dict[str, str]]:
        """Stream a headed NCBI TSV as dicts. NCBI comments out the header line ('#tax_id\\t...'), so the
        leading '#' is stripped to give a clean column name. Values are returned raw -- list splitting and
        empty-marker handling happen per-field at the call site, since which columns are lists is known."""
        logging.info(f"Streaming {path.name}..")
        with self._open_maybe_gzipped(path) as file:
            header_line = file.readline()
            if not header_line:
                return
            columns = header_line.rstrip("\n").lstrip("#").split("\t")
            for line_num, line in enumerate(file, 1):
                if line_num % 5_000_000 == 0:
                    logging.info(f"    at {line_num} rows of {path.name}")
                values = line.rstrip("\n").split("\t")
                if len(values) != len(columns):  # skip any malformed/blank line rather than mis-align columns
                    continue
                yield dict(zip(columns, values))

    def _load_taxonomy(self, input_dir: Path) -> TaxonNormalizer | None:
        """Load NCBI's taxonomy from the taxdump archive, if present. Without it, taxa stay exactly as NCBI
        Gene reports them (strain-level for many microbes) and descriptions carry no organism name."""
        taxdump_path = input_dir / TAXDUMP_FILENAME
        if taxdump_path.is_file():
            return TaxonNormalizer(taxdump_path)

        # Without the taxonomy we can't roll strain-level taxa up to species -- and the allowlist is
        # species-level, so every strain-published gene (all of yeast's, all of E. coli's) would be silently
        # dropped as "not in the allowlist". Refuse rather than quietly ingest a lopsided subset.
        if self.use_taxon_allowlist:
            raise FileNotFoundError(
                f"{self.source_name}: {TAXDUMP_FILENAME} is required in {input_dir} when use_taxon_allowlist "
                f"is on, since the allowlist is species-level and NCBI publishes many genes under strain-level "
                f"taxa. Download it from ftp.ncbi.nlm.nih.gov/pub/taxonomy/."
            )
        logging.warning(
            f"No {TAXDUMP_FILENAME} in {input_dir} - taxa will NOT be normalized to species rank and "
            f"descriptions won't name the organism. NCBI Gene reports many microbes at strain rank "
            f"(e.g. NCBITaxon:559292 rather than :4932), which won't line up with other sources. "
            f"Download it from ftp.ncbi.nlm.nih.gov/pub/taxonomy/."
        )
        return None

    def _load_summaries(self, input_dir: Path) -> dict[str, str]:
        """Build GeneID -> prose summary from gene_summary, if present. The summaries are heavily duplicated
        (3.1M rows share only ~150k distinct texts, thanks to shared RefSeq boilerplate across orthologs and
        miRNAs), so each distinct string is interned and shared -- turning ~635MB of text into ~30MB."""
        summary_path = self._resolve_optional(input_dir, SUMMARY_FILENAME)
        if not summary_path:
            logging.warning(
                f"{self.source_name}: no '{SUMMARY_FILENAME}' file in {input_dir}; genes will fall back to "
                f"NCBI's full-name 'description' column"
                if self.use_full_name_as_description_fallback
                else f"{self.source_name}: no '{SUMMARY_FILENAME}' file in {input_dir}; genes will have no description"
            )
            return {}

        summaries: dict[str, str] = {}
        interned: dict[str, str] = {}
        for row in self._stream_tsv(summary_path):
            summary = row.get(COL_SUMMARY, "").strip()
            if not summary or summary == EMPTY_VALUE:
                continue
            summaries[row[COL_GENE_ID]] = interned.setdefault(summary, summary)
        logging.info(f"Loaded {len(summaries)} gene summaries ({len(interned)} distinct texts)")
        return summaries

    # ------------------------------- row -> node -------------------------------

    def _harmonize_row(self, row: dict[str, str], summaries: dict[str, str]) -> dict[str, Any] | None:
        """Convert one gene_info row into a KRAKEN node, or None if the row is filtered out."""
        # Placeholders first: a single string compare, and doing it before the rollup keeps them out of the
        # rolled-up-to-species tally.
        if row.get(COL_SYMBOL) == PLACEHOLDER_SYMBOL:
            self.skipped_placeholders += 1
            return None

        # Then the taxon, which rejects by far the most rows. Note it tests the SPECIES-rank taxon, matching
        # how the allowlist is expressed and how `taxon` is stored.
        raw_tax_id = self._value(row, COL_TAX_ID)
        species_tax_id = self._species_tax_id(raw_tax_id)
        if self.use_taxon_allowlist and species_tax_id not in self.allowed_species:
            self.skipped_by_taxon += 1
            return None

        gene_type = self._value(row, COL_TYPE_OF_GENE)
        categories = self._categories_for(gene_type, self._values(row, COL_FEATURE_TYPE))
        if categories is None:
            self.skipped_biological_regions += 1
            return None

        gene_id = row[COL_GENE_ID]
        full_name = self._value(row, COL_AUTHORITY_FULL_NAME) or self._gene_full_name(row, gene_type)

        # Prefer the nomenclature authority's official symbol as the node name, falling back to NCBI's symbol
        name = self._value(row, COL_AUTHORITY_SYMBOL) or self._value(row, COL_SYMBOL)

        # Synonyms: every other name-ish string on the record (create_node also folds in `name` itself)
        synonyms = set(self._values(row, COL_SYNONYMS)) | set(self._values(row, COL_OTHER_DESIGNATIONS))
        synonyms |= {value for value in (full_name, self._value(row, COL_SYMBOL)) if value}

        equivalent_ids = self._equivalent_ids(row)
        summary = summaries.get(gene_id)

        # `summary or full_name` is everything we know that could describe this gene, independent of whether
        # use_full_name_as_description_fallback will actually publish the full name as the description.
        if self.drop_uncharacterized_genes and self._is_uncharacterized(
            name, self._value(row, COL_LOCUS_TAG), summary or full_name, equivalent_ids
        ):
            self.skipped_uncharacterized += 1
            return None

        description = summary
        if not description and self.use_full_name_as_description_fallback:
            description = self._full_name_description(full_name, raw_tax_id)

        taxon = f"{NCBI_TAXON_PREFIX}:{species_tax_id}" if species_tax_id else None
        # Keep the taxon exactly as NCBI reported it whenever we rolled it up, so strain-level detail survives
        # the normalization and the rollup stays auditable.
        raw_taxon_attribute = f"{NCBI_TAXON_PREFIX}:{raw_tax_id}" if species_tax_id != raw_tax_id else None

        attributes = {
            key: value
            for key, value in {
                RAW_TAXON_ATTRIBUTE: raw_taxon_attribute,
                COL_TYPE_OF_GENE: gene_type,
                COL_CHROMOSOME: self._value(row, COL_CHROMOSOME),
                COL_MAP_LOCATION: self._value(row, COL_MAP_LOCATION),
                COL_LOCUS_TAG: self._value(row, COL_LOCUS_TAG),
                COL_FEATURE_TYPE: self._values(row, COL_FEATURE_TYPE),
            }.items()
            if not is_empty(value)
        }

        return self.create_node(
            curie=f"{NCBI_GENE_PREFIX}:{gene_id}",
            categories=categories,
            provided_by=self.source_infores,
            equivalent_ids=equivalent_ids,
            name=name,
            synonyms=synonyms,
            description=description,
            urls=GENE_URL_TEMPLATE.format(gene_id=gene_id),
            taxon=taxon,
            attributes=attributes,
        )

    @staticmethod
    def _is_uncharacterized(name: str, locus_tag: str, description: str, equivalent_ids: list[str]) -> bool:
        """True only when a gene fails ALL THREE tests: it has no real symbol, nothing to cross-reference it
        to, and no description beyond a generic stand-in. Deliberately conservative -- a gene that fails just
        one or two is still reachable somehow and is kept."""
        if equivalent_ids:
            return False  # linkable, so it can still resolve against other sources
        has_real_symbol = not (LOC_PLACEHOLDER_PATTERN.match(name) or (locus_tag and name == locus_tag))
        if has_real_symbol:
            return False  # findable by name
        return not description or description.lower().startswith(GENERIC_DESCRIPTION_PREFIXES)

    @classmethod
    def _gene_full_name(cls, row: dict[str, str], gene_type: str) -> str:
        """NCBI's `description` column, which is really the gene's full name ("amyloid beta precursor
        protein"). For ~0.9% of genes it just restates type_of_gene ("pseudo", "ncRNA") -- that's a category,
        not a name, and would make a worse description than none, so it's dropped."""
        description_column = cls._value(row, COL_DESCRIPTION)
        return "" if description_column == gene_type else description_column

    def _full_name_description(self, full_name: str, tax_id: str) -> str:
        """The description used when a gene has no prose summary: its full name, qualified by the organism's
        scientific name when we know it -- e.g. "amyloid beta precursor protein (Bos taurus)". The organism
        matters here because the same gene name recurs across species, so an unqualified full name reads as
        though it were the human one. Uses the taxon exactly as NCBI reported it, so a strain-specific gene
        says so, even though the `taxon` property is rolled up to species. No taxdump loaded (or an unknown
        tax_id) just means no parenthetical."""
        taxon_name = self.taxonomy.scientific_name(tax_id) if self.taxonomy else None
        return f"{full_name} ({taxon_name})" if full_name and taxon_name else full_name

    def _categories_for(self, gene_type: str, feature_types: list[str]) -> list[str] | None:
        """Biolink categories for a record's type_of_gene. None means 'filter this record out'."""
        if gene_type in GENE_TYPES:
            return [GENE_CATEGORY]
        if gene_type == BIOLOGICAL_REGION_TYPE:
            if not self.include_biological_regions:
                return None
            is_regulatory = any(ft.startswith(REGULATORY_FEATURE_PREFIX) for ft in feature_types)
            return [REGULATORY_REGION_CATEGORY if is_regulatory else BIOLOGICAL_ENTITY_CATEGORY]
        # An unrecognized type_of_gene is still an NCBI *gene* record, so treat it as a Gene rather than
        # dropping it -- new values appear in NCBI releases and silently losing genes would be worse. Counted
        # and logged, though, so a new NCBI type doesn't quietly become a Gene without anyone noticing.
        self.unrecognized_gene_types[gene_type] += 1
        return [GENE_CATEGORY]

    def _species_tax_id(self, tax_id: str) -> str:
        """The gene's taxon rolled up to species rank. NCBI Gene reports many microbes at strain rank while
        the aggregator KGs report the species, so normalizing here keeps a `taxon` query at species level from
        silently missing strain-labeled genes. The raw strain id is kept in attributes (see _harmonize_row)."""
        if not tax_id:
            return ""
        species_tax_id = self.taxonomy.to_species(tax_id) if self.taxonomy else tax_id
        if species_tax_id != tax_id:
            self.rolled_up_to_species += 1
        return species_tax_id

    def _equivalent_ids(self, row: dict[str, str]) -> list[str]:
        """Resolve NCBI's dbXrefs into KRAKEN CURIEs via biomapper2's normalizer, which owns what counts as a
        valid prefix and id and what each one's canonical form is (so "MIM:104760" becomes OMIM:104760 and
        "Ensembl:ENSG.." becomes ENSEMBL:ENSG..). Unknown prefixes and malformed ids are counted for the
        run summary rather than silently dropped."""
        local_ids_by_prefix: dict[str, list[str]] = defaultdict(list)
        for xref in self._values(row, COL_DBXREFS):
            # NCBI writes some xrefs with the prefix doubled (e.g. "HGNC:HGNC:620", "MGI:MGI:87986")
            prefix, _, local_id = fix_repeated_prefix(xref).partition(":")
            if local_id and prefix not in DBXREF_PREFIXES_IGNORED:
                local_ids_by_prefix[prefix].append(local_id)
        if not local_ids_by_prefix:
            return []

        curies, invalid_ids, unrecognized_prefixes = self.normalizer.get_curies(
            local_ids_dict=local_ids_by_prefix,
            stop_on_invalid_id=False,
            log_warnings=False,
            fuzzy_match_vocab=False,  # exact prefix matching only; guessing at a vocab would mis-assign ids
        )
        for prefix in unrecognized_prefixes:
            self.unrecognized_dbxref_prefixes[prefix] += 1
        for prefix, bad_ids in invalid_ids.items():
            self.invalid_dbxref_ids[prefix] += len(bad_ids)
        return list(curies)

    # ------------------------------- small helpers -------------------------------

    @staticmethod
    def _value(row: dict[str, str], column: str) -> str:
        """A single scalar column value, with NCBI's '-' null marker normalized to an empty string."""
        value = row.get(column, "").strip()
        return "" if value == EMPTY_VALUE else value

    @classmethod
    def _values(cls, row: dict[str, str], column: str) -> list[str]:
        """A pipe-delimited column split into its parts (empty list for NCBI's '-' null marker)."""
        value = cls._value(row, column)
        return [part.strip() for part in value.split(LIST_DELIMITER) if part.strip()] if value else []

    def _log_run_summary(self, node_count: int):
        logging.info(f"Emitted {node_count} NCBI Gene nodes")
        if self.skipped_by_taxon:
            logging.info(
                f"Skipped {self.skipped_by_taxon} genes whose species is not among the "
                f"{len(self.allowed_species)} in ncbigene_taxon_allowlist.py"
            )
        if self.rolled_up_to_species:
            logging.info(
                f"Rolled {self.rolled_up_to_species} genes' taxa up to species rank (NCBI reported them at "
                f"strain/subspecies rank); the reported taxon is kept in the '{RAW_TAXON_ATTRIBUTE}' attribute"
            )
        if self.skipped_uncharacterized:
            logging.info(
                f"Skipped {self.skipped_uncharacterized} uncharacterized genes (no real symbol, no dbXref, and "
                f"no description beyond a generic stand-in); set drop_uncharacterized_genes=False to keep them"
            )
        if self.skipped_placeholders:
            logging.info(f"Skipped {self.skipped_placeholders} '{PLACEHOLDER_SYMBOL}' placeholder records")
        if self.skipped_biological_regions:
            logging.info(
                f"Skipped {self.skipped_biological_regions} '{BIOLOGICAL_REGION_TYPE}' records "
                f"(set include_biological_regions=True to ingest them)"
            )
        if self.unrecognized_gene_types:
            ranked = dict(sorted(self.unrecognized_gene_types.items(), key=lambda kv: kv[1], reverse=True))
            logging.warning(
                f"Typed as {GENE_CATEGORY} by fallback: {len(ranked)} type_of_gene value(s) not in GENE_TYPES "
                f"- counts are: {ranked}. Check whether {GENE_CATEGORY} is right for these, and add them to "
                f"GENE_TYPES (or give them their own category) either way."
            )
        if self.unrecognized_dbxref_prefixes:
            ranked = dict(sorted(self.unrecognized_dbxref_prefixes.items(), key=lambda kv: kv[1], reverse=True))
            logging.warning(
                f"Dropped dbXrefs from {len(ranked)} prefix(es) biomapper2 doesn't recognize - counts by "
                f"prefix are: {ranked}. Register any that are 1:1 with an NCBI gene in biomapper2 "
                f"(core/normalizer/vocab_config.py), or add it to DBXREF_PREFIXES_IGNORED here if it isn't a "
                f"gene identity space."
            )
        if self.invalid_dbxref_ids:
            ranked = dict(sorted(self.invalid_dbxref_ids.items(), key=lambda kv: kv[1], reverse=True))
            logging.warning(f"Dropped dbXrefs whose ids failed biomapper2 validation - counts by prefix are: {ranked}")
