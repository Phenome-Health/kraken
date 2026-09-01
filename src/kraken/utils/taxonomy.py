"""NCBI taxonomy lookups, backed by NCBI's taxdump archive (ftp.ncbi.nlm.nih.gov/pub/taxonomy/).

Two jobs, both harmonization-time concerns:
  - naming a taxon (`scientific_name`), so a gene's description can say which organism it belongs to
  - rolling a taxon up to species rank (`to_species`)

The rollup exists because sources cite the taxonomy at different DEPTHS. NCBI Gene labels many microbes by
strain -- S. cerevisiae genes carry NCBITaxon:559292 (S288C), E. coli genes NCBITaxon:511145 (K-12 MG1655) --
while the aggregator KGs cite the species (4932, 562). Those are the same organism written two ways, so
without normalization a `taxon` query at species level silently misses the strain-labeled nodes, and any
one-taxon-per-entity merge rule would treat them as a conflict and refuse a correct merge.

Normalization only ever moves UP the tree. A rank below species (strain, subspecies, serovar, ...) becomes its
species ancestor; a species stays put; anything at or above species, or with no species ancestor at all
(many viruses, "unclassified" clades), is returned unchanged -- there is no defensible way to push a genus
down to a species. Strain detail isn't lost: callers are expected to retain the raw taxon alongside.
"""

import logging
import tarfile
from pathlib import Path

NODES_DMP_MEMBER = "nodes.dmp"
NAMES_DMP_MEMBER = "names.dmp"

# Both .dmp files are 'field\t|\tfield\t|\t...' rows
DMP_DELIMITER = "\t|"

NODES_TAX_ID_FIELD = 0
NODES_PARENT_FIELD = 1
NODES_RANK_FIELD = 2

NAMES_TAX_ID_FIELD = 0
NAMES_NAME_FIELD = 1
NAMES_CLASS_FIELD = 3
SCIENTIFIC_NAME_CLASS = "scientific name"

SPECIES_RANK = "species"
ROOT_TAX_ID = "1"

NCBI_TAXON_PREFIX = "NCBITaxon"


class TaxonNormalizer:
    """Rolls taxa up to species rank and looks up their scientific names.

    Built from a taxdump archive, read in place (never extracted to disk). Only the taxa that actually need
    rolling up are retained -- most of NCBI's ~2.7M taxa are already species rank and map to themselves, so
    the resident map holds just the sub-species entries.
    """

    def __init__(self, taxdump_path: Path, include_names: bool = True):
        self.taxdump_path = Path(taxdump_path)
        self._species_of: dict[str, str] = {}  # only sub-species taxa; everything else maps to itself
        self._names: dict[str, str] = {}
        self._load(include_names)

    # ------------------------------- public API -------------------------------

    def to_species(self, tax_id: str) -> str:
        """The species-rank ancestor of `tax_id`, or `tax_id` itself when it is a species, is above species,
        or has no species ancestor. Accepts a bare id ('559292') or a CURIE ('NCBITaxon:559292') and returns
        the same form it was given."""
        prefix, bare_id = self._split(tax_id)
        species = self._species_of.get(bare_id, bare_id)
        return f"{prefix}:{species}" if prefix else species

    def scientific_name(self, tax_id: str) -> str | None:
        """The taxon's scientific name, or None if unknown (or names weren't loaded)."""
        return self._names.get(self._split(tax_id)[1])

    # ------------------------------- loading -------------------------------

    @staticmethod
    def _split(tax_id: str) -> tuple[str | None, str]:
        """Split an optional 'NCBITaxon:' prefix off, so callers can pass either form."""
        prefix, delimiter, bare_id = tax_id.partition(":")
        return (prefix, bare_id) if delimiter else (None, tax_id)

    def _load(self, include_names: bool) -> None:
        if not self.taxdump_path.is_file():
            raise FileNotFoundError(f"No taxdump archive at {self.taxdump_path}")
        with tarfile.open(self.taxdump_path, "r:gz") as archive:
            parents, species_taxa = self._parse_nodes(archive)
            if include_names:
                self._names = self._parse_names(archive)
        self._species_of = self._build_species_map(parents, species_taxa)
        logging.info(
            f"Loaded taxonomy from {self.taxdump_path.name}: {len(parents)} taxa, "
            f"{len(self._species_of)} of them below species rank"
        )

    def _parse_nodes(self, archive: tarfile.TarFile) -> tuple[dict[str, str], set[str]]:
        """Read nodes.dmp into (tax_id -> parent tax_id, set of species-rank tax_ids)."""
        node_file = archive.extractfile(NODES_DMP_MEMBER)
        if node_file is None:
            raise ValueError(f"{self.taxdump_path} has no '{NODES_DMP_MEMBER}' member")
        parents: dict[str, str] = {}
        species_taxa: set[str] = set()
        for raw_line in node_file:
            fields = [field.strip() for field in raw_line.decode("utf-8").split(DMP_DELIMITER)]
            if len(fields) <= NODES_RANK_FIELD:
                continue
            tax_id = fields[NODES_TAX_ID_FIELD]
            parents[tax_id] = fields[NODES_PARENT_FIELD]
            if fields[NODES_RANK_FIELD] == SPECIES_RANK:
                species_taxa.add(tax_id)
        return parents, species_taxa

    def _parse_names(self, archive: tarfile.TarFile) -> dict[str, str]:
        """Read names.dmp into tax_id -> scientific name (each taxon also carries synonyms and common names,
        which we skip)."""
        names_file = archive.extractfile(NAMES_DMP_MEMBER)
        if names_file is None:
            raise ValueError(f"{self.taxdump_path} has no '{NAMES_DMP_MEMBER}' member")
        names: dict[str, str] = {}
        for raw_line in names_file:
            fields = [field.strip() for field in raw_line.decode("utf-8").split(DMP_DELIMITER)]
            if len(fields) > NAMES_CLASS_FIELD and fields[NAMES_CLASS_FIELD] == SCIENTIFIC_NAME_CLASS:
                names[fields[NAMES_TAX_ID_FIELD]] = fields[NAMES_NAME_FIELD]
        return names

    @staticmethod
    def _build_species_map(parents: dict[str, str], species_taxa: set[str]) -> dict[str, str]:
        """Precompute sub-species tax_id -> species ancestor, walking each taxon up until a species-rank
        ancestor is found (or the walk leaves the tree). Only taxa that actually resolve to a DIFFERENT id are
        kept, so species and above -- the large majority -- cost nothing.

        Walks are memoized against the map being built, so shared lineages are traversed once.
        """
        species_of: dict[str, str] = {}
        for start_tax_id in parents:
            if start_tax_id in species_taxa:
                continue  # already a species; maps to itself
            path: list[str] = []
            current = start_tax_id
            # Climb until we hit a species, a taxon we've already resolved, or the top of the tree. The
            # `current in path` guard is belt-and-braces against a malformed dump introducing a cycle.
            while (
                current not in species_taxa
                and current not in species_of
                and current != ROOT_TAX_ID
                and current in parents
                and current not in path
            ):
                path.append(current)
                current = parents[current]

            resolved = current if current in species_taxa else species_of.get(current)
            if resolved:
                for tax_id in path:
                    species_of[tax_id] = resolved
        return species_of
