# bio_bmi.py
import csv
import logging
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import (
    BIOLOGICAL_BMI_SOURCE_ID,
    DATA_ANALYSIS_PIPELINE,
    KNOWLEDGE_ASSERTION,
    MANUAL_AGENT,
    OBJ_DIRECTION_QUALIFIER,
    STATISTICAL_ASSOCIATION,
)
from kraken.utils.kg_io import save_to_jsonl

PAPER = "PMID:36941332"  # Watanabe et al., Nat Med 2023 -- doi:10.1038/s41591-023-02248-0

# --- PLACEHOLDER Biolink types / minted ids (pending a types review; each a single swappable constant) ---
# A "biological BMI" score is a measurable characteristic of a person, so Attribute (not InformationContentEntity).
NODE_CATEGORY = "biolink:Attribute"
COMPONENT_STUB_CATEGORY = "biolink:NamedThing"  # real type arrives when the CHEBI/NCBIGene/etc node merges in
MODEL_OF_PREDICATE = "biolink:model_of"  # model -> Biological BMI (the model approximates it)
CONTRIBUTES_TO_PREDICATE = "biolink:contributes_to"  # component -> model (the feature contributes to the score)
BIO_BMI_PREFIX = "BIOBMI"  # minted ids for the umbrella + model nodes (no registered ontology for these)

UMBRELLA_ID = f"{BIO_BMI_PREFIX}:BiologicalBMI"

# model id -> (display name, description, TSV filename)
MODELS = {
    f"{BIO_BMI_PREFIX}:MetBMI": (
        "Metabolomic BMI",
        f"Metabolomics-based model estimating biological BMI (Watanabe et al., Nat Med 2023).",
        ["MetBMI"],
        "metbmi_entity_node_map_SHARED.tsv",
    ),
    f"{BIO_BMI_PREFIX}:ProtBMI": (
        "Proteomic BMI",
        f"Proteomics-based model estimating biological BMI (Watanabe et al., Nat Med 2023).",
        ["ProtBMI"],
        "protbmi_entity_node_map_SHARED.tsv",
    ),
    f"{BIO_BMI_PREFIX}:CombiBMI": (
        "Combined multiomic BMI",
        f"Combined multiomic (metabolomics, proteomics, clinical labs) model estimating biological "
        f"BMI (Watanabe et al., Nat Med 2023).",
        ["CombiBMI"],
        "combi_entity_node_map_SHARED.tsv",
    ),
}


class BioBMIHarmonizer(BaseHarmonizer):
    """Harmonizer for the multiomic 'biological BMI' models (Watanabe et al., Nat Med 2023;
    doi:10.1038/s41591-023-02248-0). This is a paper-derived source, not a database.

    Builds an umbrella 'Biological BMI' node plus three model nodes (MetBMI/ProtBMI/CombiBMI), each linked to the
    umbrella via `model_of`, and links every model to its component omics features (metabolites / proteins /
    clinical labs) via `contributes_to`. Components come from per-model TSVs whose entities were pre-mapped to
    kraken CURIEs; we use the primary `kraken_id` only (the alt2-5 columns are unreliable -- for proteins they
    hold *other* genes), normalized via biomapper2, and stream nodes/edges per model (components deduped across
    models). Each feature's direction of association with BMI is an object_direction_qualifier; the model
    coefficient (mean_beta) and omics_type are edge attributes.

    Node/edge Biolink types are placeholders -- see the module-level constants; each is a single swappable value.
    """

    source_infores = BIOLOGICAL_BMI_SOURCE_ID

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self._curie_cache: dict[str, str | None] = {}  # raw kraken_id -> normalized curie (or None)

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
        input_dir = Path(input_file)
        logging.info(f"Harmonizing {self.source_name}: {input_dir} -> {nodes_output}, {edges_output}")

        save_to_jsonl([], nodes_output, mode="w")  # truncate; we stream per-model below
        save_to_jsonl([], edges_output, mode="w")

        # 1. Umbrella + model nodes, and model_of edges (model -> Biological BMI)
        top_nodes = [
            self._make_attr_node(
                UMBRELLA_ID,
                "Biological BMI",
                f"Multiomic estimate of body mass index reflecting metabolic health.",
                ["BioBMI", "bBMI"],
            )
        ]
        top_edges = []
        for model_curie, (name, description, synonyms, _file) in MODELS.items():
            top_nodes.append(self._make_attr_node(model_curie, name, description, synonyms))
            top_edges.append(self._make_model_of_edge(model_curie))
        save_to_jsonl(top_nodes, nodes_output, mode="a")
        save_to_jsonl(top_edges, edges_output, mode="a")

        # 2. Per-model component nodes + contributes_to edges (streamed; components deduped across models)
        written_components: set[str] = set()
        total_edges = total_unmapped = 0
        for model_curie, (name, _desc, _syn, filename) in MODELS.items():
            node_batch, edge_batch = [], []
            for row in self._read_tsv(input_dir / filename):
                comp_curie = self._normalize_component(row.get("kraken_id"))
                if not comp_curie:
                    total_unmapped += 1
                    continue
                if comp_curie not in written_components:
                    written_components.add(comp_curie)
                    node_batch.append(self._make_component_stub(comp_curie, row.get("kraken_name")))
                edge_batch.append(self._make_contributes_edge(comp_curie, model_curie, row))
            save_to_jsonl(node_batch, nodes_output, mode="a")
            save_to_jsonl(edge_batch, edges_output, mode="a")
            total_edges += len(edge_batch)
            logging.info(f"  {name}: {len(edge_batch)} component edges")

        logging.info(
            f"{self.source_name} harmonization complete: {1 + len(MODELS)} model nodes + "
            f"{len(written_components)} component stubs; {len(MODELS)} model_of + {total_edges} contributes_to "
            f"edges ({total_unmapped} rows skipped -- unmapped/unresolvable)"
        )

    # ------------------------------------------------------------------ readers

    @staticmethod
    def _read_tsv(path: Path):
        with open(path, encoding="utf-8") as f:
            yield from csv.DictReader(f, delimiter="\t")

    # ------------------------------------------------------------------ node/edge builders

    def _make_attr_node(self, curie: str, name: str, description: str, synonyms: list[str]) -> dict:
        return self.create_node(
            curie=curie,
            categories=[NODE_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=name,
            synonyms=[name] + synonyms,
            description=description,
            publications=[PAPER],
        )

    def _make_component_stub(self, curie: str, name: str | None) -> dict:
        return self.create_node(
            curie=curie,
            categories=[COMPONENT_STUB_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=(name or None),
        )

    def _make_model_of_edge(self, model_curie: str) -> dict:
        return self.create_edge(
            subject_id=model_curie,
            object_id=UMBRELLA_ID,
            predicate=MODEL_OF_PREDICATE,
            primary_ks=self.source_infores,
            knowledge_level=KNOWLEDGE_ASSERTION,
            agent_type=MANUAL_AGENT,
            publications=[PAPER],
        )

    def _make_contributes_edge(self, comp_curie: str, model_curie: str, row: dict) -> dict:
        direction = (row.get("direction") or "").strip()
        qualifiers = {OBJ_DIRECTION_QUALIFIER: direction} if direction in ("increased", "decreased") else None
        attributes = {"omics_type": (row.get("omics_type") or "").strip() or None}
        mean_beta = self._to_float(row.get("mean_beta"))
        if mean_beta is not None:
            attributes["mean_beta"] = round(mean_beta, 6)
        return self.create_edge(
            subject_id=comp_curie,
            object_id=model_curie,
            predicate=CONTRIBUTES_TO_PREDICATE,
            primary_ks=self.source_infores,
            knowledge_level=STATISTICAL_ASSOCIATION,
            agent_type=DATA_ANALYSIS_PIPELINE,
            publications=[PAPER],
            qualifiers=qualifiers,
            attributes={k: v for k, v in attributes.items() if v not in (None, "")},
        )

    # ------------------------------------------------------------------ helpers

    def _normalize_component(self, raw_id: str | None) -> str | None:
        """Validate/canonicalize a pre-mapped kraken_id (a full CURIE like 'CHEBI:64583') via biomapper2.
        Returns None (and warns once) for empty ids or vocabs biomapper doesn't recognize (e.g. PathWhiz)."""
        raw_id = (raw_id or "").strip()
        if not raw_id or ":" not in raw_id:
            return None
        if raw_id in self._curie_cache:
            return self._curie_cache[raw_id]
        prefix, _, local = raw_id.partition(":")
        resolved, _, _ = self.normalizer.get_curies({prefix: local}, stop_on_invalid_id=False, log_warnings=False)
        curie = next(iter(resolved), None)
        if not curie:
            logging.warning(f"{self.source_name}: could not normalize component id {raw_id!r}; dropping it.")
        self._curie_cache[raw_id] = curie
        return curie

    @staticmethod
    def _to_float(value) -> float | None:
        try:
            return float(value) if value not in (None, "") else None
        except (ValueError, TypeError):
            return None
