# bio_age.py
import csv
import logging
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import (
    BIO_AGE_SOURCE_ID,
    DATA_ANALYSIS_PIPELINE,
    STATISTICAL_ASSOCIATION,
)
from kraken.utils.kg_io import save_to_jsonl

PAPER = "PMID:31724055"  # Earls et al., J Gerontol A 2019 -- doi:10.1093/gerona/glz220

# --- Biolink types / minted ids (each a single swappable constant) ---
# A biological-age score is a measurable characteristic of a person, so Attribute.
NODE_CATEGORY = "biolink:Attribute"
CORRELATE_STUB_CATEGORY = "biolink:NamedThing"  # disease/behavior endpoints; real type arrives on merge
# ΔAge direction selects the predicate (correlation direction IS the predicate, so no qualifier needed):
POSITIVELY_CORRELATED = "biolink:positively_correlated_with"  # condition -> increased ΔAge
NEGATIVELY_CORRELATED = "biolink:negatively_correlated_with"  # condition -> decreased ΔAge
BIO_AGE_PREFIX = "BIOAGE"  # minted id for the Biological Age node (no registered ontology for it)

BIO_AGE_ID = f"{BIO_AGE_PREFIX}:BiologicalAge"

# Curated TSV columns. Only kraken_curie + delta_age_direction are required; the ΔAge magnitude/CI columns are
# optional and become edge attributes when present (pending the numeric Fig 2 values from the paper's authors).
COL_CURIE = "kraken_curie"
COL_DIRECTION = "delta_age_direction"
COL_NAME = "kraken_name"
COL_DELTA_AGE = "delta_age_years"
COL_CI_LOW = "ci_low"
COL_CI_HIGH = "ci_high"

MAPPING_FILENAME = "bioage_condition_mappings.tsv"


class BioAgeHarmonizer(BaseHarmonizer):
    """Harmonizer for multi-omic biological age (bA) (Earls et al., J Gerontol A 2019;
    doi:10.1093/gerona/glz220). A paper-derived source, not a database.

    Builds a single 'Biological Age' node and links it to the health conditions/behaviors it correlates with,
    taken from Figure 2 of the paper (ΔAge associations). We deliberately do NOT model bA's input analytes (the
    Klemera-Doubal models are PCA-based, so *every* analyte gets a nonzero coefficient -- there is no honest
    sparsity to model). Each condition was pre-mapped to a kraken CURIE (curated `kraken_curie` column); the ΔAge
    direction selects the predicate (increased -> positively_correlated_with, decreased ->
    negatively_correlated_with). ΔAge magnitude + CI columns, if present, become edge attributes.

    Node/edge Biolink types are single swappable module-level constants -- see above.
    """

    source_infores = BIO_AGE_SOURCE_ID

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self._curie_cache: dict[str, str | None] = {}  # raw curie -> normalized curie (or None)

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
            raise ValueError(f"{self.source_name} requires input_file (the bio-age directory)")
        mapping_path = Path(input_file) / MAPPING_FILENAME
        logging.info(f"Harmonizing {self.source_name}: {mapping_path} -> {nodes_output}, {edges_output}")

        save_to_jsonl([], nodes_output, mode="w")  # truncate; we stream below
        save_to_jsonl([], edges_output, mode="w")

        # 1. The single Biological Age node
        save_to_jsonl(
            [
                self._make_attr_node(
                    BIO_AGE_ID,
                    "Biological Age",
                    "Multi-omic Klemera-Doubal estimate of biological age (Earls et al., J Gerontol A 2019).",
                )
            ],
            nodes_output,
            mode="a",
        )

        # 2. Condition correlation edges (Biological Age -> condition), streamed; endpoints deduped
        written_stubs: set[str] = set()
        node_batch, edge_batch = [], []
        total_unmapped = skipped_dir = 0
        for row in self._read_tsv(mapping_path):
            predicate = self._predicate_for_direction(row.get(COL_DIRECTION))
            if predicate is None:
                skipped_dir += 1
                continue
            curie = self._normalize_curie(row.get(COL_CURIE))
            if not curie:
                total_unmapped += 1
                continue
            if curie not in written_stubs:
                written_stubs.add(curie)
                node_batch.append(self._make_stub(curie, row.get(COL_NAME)))
            edge_batch.append(self._make_correlation_edge(curie, predicate, row))
        save_to_jsonl(node_batch, nodes_output, mode="a")
        save_to_jsonl(edge_batch, edges_output, mode="a")

        logging.info(
            f"{self.source_name} harmonization complete: 1 Biological Age node + {len(written_stubs)} condition "
            f"stubs; {len(edge_batch)} correlation edges "
            f"({total_unmapped} unmapped, {skipped_dir} skipped for missing/invalid direction)"
        )

    # ------------------------------------------------------------------ readers

    @staticmethod
    def _read_tsv(path: Path):
        with open(path, encoding="utf-8") as f:
            yield from csv.DictReader(f, delimiter="\t")

    # ------------------------------------------------------------------ node/edge builders

    def _make_attr_node(self, curie: str, name: str, description: str) -> dict:
        return self.create_node(
            curie=curie,
            categories=[NODE_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=name,
            description=description,
            publications=[PAPER],
        )

    def _make_stub(self, curie: str, name: str | None) -> dict:
        return self.create_node(
            curie=curie,
            categories=[CORRELATE_STUB_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=(name or None),
        )

    def _make_correlation_edge(self, curie: str, predicate: str, row: dict) -> dict:
        attributes = {}
        for col in (COL_DELTA_AGE, COL_CI_LOW, COL_CI_HIGH):
            value = self._to_float(row.get(col))
            if value is not None:
                attributes[col] = value
        return self.create_edge(
            subject_id=BIO_AGE_ID,
            object_id=curie,
            predicate=predicate,
            primary_ks=self.source_infores,
            knowledge_level=STATISTICAL_ASSOCIATION,
            agent_type=DATA_ANALYSIS_PIPELINE,
            publications=[PAPER],
            attributes=attributes or None,
        )

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _predicate_for_direction(direction: str | None) -> str | None:
        d = (direction or "").strip().lower()
        if d == "increased":
            return POSITIVELY_CORRELATED
        if d == "decreased":
            return NEGATIVELY_CORRELATED
        return None

    def _normalize_curie(self, raw_id: str | None) -> str | None:
        """Validate/canonicalize a curated kraken CURIE via biomapper2. Returns None (and warns once) for empty
        ids or vocabs biomapper2 doesn't recognize."""
        raw_id = (raw_id or "").strip()
        if not raw_id or ":" not in raw_id:
            return None
        if raw_id in self._curie_cache:
            return self._curie_cache[raw_id]
        prefix, _, local = raw_id.partition(":")
        resolved, _, _ = self.normalizer.get_curies({prefix: local}, stop_on_invalid_id=False, log_warnings=False)
        curie = next(iter(resolved), None)
        if not curie:
            logging.warning(f"{self.source_name}: could not normalize {raw_id!r}; dropping it.")
        self._curie_cache[raw_id] = curie
        return curie

    @staticmethod
    def _to_float(value) -> float | None:
        try:
            return float(value) if value not in (None, "") else None
        except (ValueError, TypeError):
            return None
