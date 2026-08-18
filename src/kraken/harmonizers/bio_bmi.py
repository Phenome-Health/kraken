# bio_bmi.py
import csv
import logging
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import (
    BIOLOGICAL_BMI_SOURCE_ID,
    DATA_ANALYSIS_PIPELINE,
    STATISTICAL_ASSOCIATION,
)
from kraken.utils.kg_io import save_to_jsonl

PAPER = "PMID:36941332"  # Watanabe et al., Nat Med 2023 -- doi:10.1038/s41591-023-02248-0

# --- Biolink types / minted ids (each a single swappable constant) ---
# A biological-BMI score is a measurable characteristic of a person, so Attribute.
NODE_CATEGORY = "biolink:Attribute"
CORRELATE_STUB_CATEGORY = "biolink:NamedThing"  # physiological-feature endpoints; real type arrives on merge
# Coefficient sign selects the predicate (correlation direction IS the predicate, so no qualifier needed):
POSITIVELY_CORRELATED = "biolink:positively_correlated_with"  # feature -> higher biological BMI
NEGATIVELY_CORRELATED = "biolink:negatively_correlated_with"  # feature -> lower biological BMI
BIO_BMI_PREFIX = "BIOBMI"  # minted id for the Biological BMI node (no registered ontology for it)

BIO_BMI_ID = f"{BIO_BMI_PREFIX}:BiologicalBMI"

# Curated TSV: physiological features regressed against several BMI estimates, with a manually-added
# KRAKEN_node_id mapping column. We keep only rows that were mapped AND are significantly correlated with the
# combined multiomic biological BMI (CombiBMI). Direction comes from the CombiBMI coefficient's sign; every
# model's coefficient / 95% CI / adjusted p-value is stashed on the edge as attributes.
MAPPING_FILENAME = "PhysiologicalFeatures-Table 1.tsv"
COL_CURIE = "KRAKEN_node_id"
COL_NAME = "KRAKEN_node_name"

# The BMI estimates present in the table, in column-prefix form. "BMI" is actual measured BMI; the rest are the
# omics/clinical models from Watanabe et al. CombiBMI (combined multiomic) is treated as the overall biological BMI.
MODEL_COL_PREFIXES = ["BMI", "MetBMI", "ProtBMI", "ChemBMI", "CombiBMI"]
PRIMARY_MODEL = "CombiBMI"  # drives edge direction + the significance filter
SIG_THRESHOLD = 0.05  # keep rows with {PRIMARY_MODEL}_AdjPval_all below this


class BioBMIHarmonizer(BaseHarmonizer):
    """Harmonizer for multiomic 'biological BMI' (Watanabe et al., Nat Med 2023;
    doi:10.1038/s41591-023-02248-0). A paper-derived source, not a database.

    Builds a single 'Biological BMI' node and links it to the physiological features it correlates with. We
    deliberately do NOT model biological BMI's input analytes (the models are regularized/omics-wide, so there is
    no honest sparsity to model). Instead we use a curated table of physiological features (waist-to-height ratio,
    blood pressure, activity metrics, polygenic scores, etc.) each regressed against the BMI estimates; rows carry
    a manually-added KRAKEN_node_id. We keep rows that are mapped AND significantly correlated with the combined
    multiomic biological BMI (CombiBMI adjusted p < 0.05). The CombiBMI coefficient's sign selects the predicate
    (positive -> positively_correlated_with, negative -> negatively_correlated_with); every model's coefficient,
    95% CI, and adjusted p-value are attached as edge attributes.

    Node/edge Biolink types are single swappable module-level constants -- see above.
    """

    source_infores = BIOLOGICAL_BMI_SOURCE_ID

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
            raise ValueError(f"{self.source_name} requires input_file (the bio-bmi directory)")
        mapping_path = Path(input_file) / MAPPING_FILENAME
        logging.info(f"Harmonizing {self.source_name}: {mapping_path} -> {nodes_output}, {edges_output}")

        save_to_jsonl([], nodes_output, mode="w")  # truncate; we stream below
        save_to_jsonl([], edges_output, mode="w")

        # 1. The single Biological BMI node
        save_to_jsonl(
            [
                self._make_attr_node(
                    BIO_BMI_ID,
                    "Biological BMI",
                    "Multiomic estimate of body mass index reflecting metabolic health "
                    "(Watanabe et al., Nat Med 2023).",
                )
            ],
            nodes_output,
            mode="a",
        )

        # 2. Feature correlation edges (Biological BMI -> feature), streamed; endpoints deduped
        written_stubs: set[str] = set()
        node_batch, edge_batch = [], []
        n_unmapped = n_nonsig = n_no_direction = 0
        for row in self._read_tsv(mapping_path):
            if not (row.get(COL_CURIE) or "").strip():
                n_unmapped += 1
                continue
            if not self._is_significant(row):
                n_nonsig += 1
                continue
            predicate = self._predicate_for_coef(row.get(f"{PRIMARY_MODEL}_Bcoef"))
            if predicate is None:
                n_no_direction += 1
                continue
            curie = self._normalize_curie(row.get(COL_CURIE))
            if not curie:
                n_unmapped += 1
                continue
            if curie not in written_stubs:
                written_stubs.add(curie)
                node_batch.append(self._make_stub(curie, row.get(COL_NAME)))
            edge_batch.append(self._make_correlation_edge(curie, predicate, row))
        save_to_jsonl(node_batch, nodes_output, mode="a")
        save_to_jsonl(edge_batch, edges_output, mode="a")

        logging.info(
            f"{self.source_name} harmonization complete: 1 Biological BMI node + {len(written_stubs)} feature "
            f"stubs; {len(edge_batch)} correlation edges ({n_unmapped} unmapped/unresolvable, {n_nonsig} "
            f"non-significant, {n_no_direction} missing coefficient)"
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
        attributes: dict[str, float] = {}
        for prefix in MODEL_COL_PREFIXES:
            key = prefix.lower()
            bcoef = self._to_float(row.get(f"{prefix}_Bcoef"))
            if bcoef is not None:
                attributes[f"{key}_bcoef"] = round(bcoef, 6)
            ci_low = self._to_float(row.get(f"{prefix}_BcoefCIlow"))
            if ci_low is not None:
                attributes[f"{key}_bcoef_ci_low"] = round(ci_low, 6)
            ci_high = self._to_float(row.get(f"{prefix}_BcoefCIhigh"))
            if ci_high is not None:
                attributes[f"{key}_bcoef_ci_high"] = round(ci_high, 6)
            adj_pval = self._to_float(row.get(f"{prefix}_AdjPval_all"))
            if adj_pval is not None:
                attributes[f"{key}_adj_pval"] = adj_pval  # not rounded: p-values can be extremely small
        return self.create_edge(
            subject_id=BIO_BMI_ID,
            object_id=curie,
            predicate=predicate,
            primary_ks=self.source_infores,
            knowledge_level=STATISTICAL_ASSOCIATION,
            agent_type=DATA_ANALYSIS_PIPELINE,
            publications=[PAPER],
            attributes=attributes or None,
        )

    # ------------------------------------------------------------------ helpers

    def _is_significant(self, row: dict) -> bool:
        pval = self._to_float(row.get(f"{PRIMARY_MODEL}_AdjPval_all"))
        return pval is not None and pval < SIG_THRESHOLD

    @staticmethod
    def _predicate_for_coef(bcoef: str | None) -> str | None:
        value = BioBMIHarmonizer._to_float(bcoef)
        if value is None or value == 0:
            return None
        return POSITIVELY_CORRELATED if value > 0 else NEGATIVELY_CORRELATED

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
