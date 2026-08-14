# loinc.py
import csv
import logging
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import LOINC_INFORES
from kraken.utils.kg_io import save_to_jsonl

# Quick-and-dirty LOINC terminology ingest: one node per LOINC identifier, with a display name + name synonyms.
# Nodes only (no edges). Everything lives under the single LOINC: prefix (numeric codes, LP parts, LA answers,
# LL answer lists, LG groups).
LOINC_PREFIX = "LOINC"

# Per-identifier-type Biolink categories. LOINC is the source of truth for most of these, so we assign
# meaningful types up front rather than defaulting to NamedThing.
#   - Numeric codes are heterogeneous: genuine measurements (quantitative/ordinal/nominal lab & clinical
#     results) -> ClinicalMeasurement; documents, narratives, and survey/claims items -> InformationContentEntity;
#     panels/sets & anything else -> NamedThing. Read off SCALE_TYP + CLASSTYPE per row (see numeric_category).
#   - LA answers / LL answer lists / LG groups are information artifacts -> InformationContentEntity.
#   - LP parts are genuinely heterogeneous (components, specimens, methods, genes, radiology axes...) with no
#     single honest fit, so they stay NamedThing -- leaf-filtering upgrades the ones that merge with a richer
#     node (e.g. the 'Creatinine' component part merges into the CHEBI ChemicalEntity node).
MEASUREMENT_CATEGORY = "biolink:ClinicalMeasurement"
INFO_CATEGORY = "biolink:InformationContentEntity"
NAMED_THING = "biolink:NamedThing"

MEASUREMENT_SCALES = {"Qn", "SemiQn", "OrdQn", "Nom", "Ord"}  # LOINC SCALE_TYP values denoting a measurement
NARRATIVE_SCALES = {"Doc", "Nar"}  # documents / narrative reports
INFO_CLASSTYPES = {"3", "4"}  # LOINC CLASSTYPE 3=Claims/Attachments, 4=Survey (questionnaire items)


def numeric_category(row: dict) -> str:
    """Category for a numeric LOINC code, from its SCALE_TYP + CLASSTYPE (see the mapping notes above)."""
    if (row.get("CLASSTYPE") or "").strip() in INFO_CLASSTYPES:
        return INFO_CATEGORY
    scale = (row.get("SCALE_TYP") or "").strip()
    if scale in NARRATIVE_SCALES:
        return INFO_CATEGORY
    if scale in MEASUREMENT_SCALES:
        return MEASUREMENT_CATEGORY
    return NAMED_THING


# "Simple" files: one identifier per row. (relative path, id col, name col, [synonym cols], category).
# `category` is a fixed string, or a callable(row) -> category for per-row typing (numeric codes).
SIMPLE_FILES = [
    ("LoincTable/Loinc.csv", "LOINC_NUM", "LONG_COMMON_NAME", ["SHORTNAME", "DisplayName", "CONSUMER_NAME"], numeric_category),
    ("AccessoryFiles/PartFile/Part.csv", "PartNumber", "PartName", ["PartDisplayName"], NAMED_THING),
    ("AccessoryFiles/GroupFile/Group.csv", "GroupId", "Group", [], INFO_CATEGORY),
    ("AccessoryFiles/GroupFile/ParentGroup.csv", "ParentGroupId", "ParentGroup", [], INFO_CATEGORY),
]

# Answer file: each row carries two identifiers -- an answer (LA...) and its answer list (LL...). (id, name, category)
ANSWER_FILE = "AccessoryFiles/PanelsAndForms/AnswerList.csv"
ANSWER_ID_NAME_COLS = [("AnswerStringId", "DisplayText", INFO_CATEGORY), ("AnswerListId", "AnswerListName", INFO_CATEGORY)]

BATCH_SIZE = 50_000


class LoincHarmonizer(BaseHarmonizer):
    """Quick-and-dirty LOINC ingest. Mints one node per LOINC identifier (numeric lab codes, LP parts, LA
    answers, LL answer lists, LG groups) with a display name and any short/display-name synonyms, all under the
    LOINC: prefix. Nodes only -- no edges, no hierarchy. Node category is a single swappable placeholder
    (NamedThing); the real type is picked up when the corresponding node from a richer source merges in.
    """

    source_infores = LOINC_INFORES

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self._curie_cache: dict[str, str | None] = {}  # LOINC local id -> canonical curie (or None)
        self._n_unmapped = 0

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
            raise ValueError(f"{self.source_name} requires input_file (the LOINC release directory)")
        base = Path(input_file)
        logging.info(f"Harmonizing {self.source_name}: {base} -> {nodes_output}")

        save_to_jsonl([], nodes_output, mode="w")  # truncate; we stream below
        save_to_jsonl([], edges_output, mode="w")  # nodes-only source

        self._n_unmapped = 0
        seen: set[str] = set()
        counts: dict[str, int] = {}
        for rel_path, id_col, name_col, syn_cols, category in SIMPLE_FILES:
            counts[rel_path] = self._ingest(base / rel_path, [(id_col, name_col, category)], syn_cols, nodes_output, seen)
        counts[ANSWER_FILE] = self._ingest(base / ANSWER_FILE, ANSWER_ID_NAME_COLS, [], nodes_output, seen)

        summary = ", ".join(f"{Path(p).name}={n}" for p, n in counts.items())
        logging.info(
            f"{self.source_name} harmonization complete: {len(seen)} LOINC nodes ({summary})"
            + (f"; {self._n_unmapped} ids dropped (biomapper2 could not form a CURIE)" if self._n_unmapped else "")
        )

    # ------------------------------------------------------------------ core

    def _ingest(
        self,
        path: Path,
        id_name_cats: list[tuple[str, str, str]],
        syn_cols: list[str],
        nodes_output: Path,
        seen: set[str],
    ) -> int:
        """Stream a LOINC CSV, minting a node for each (id_col, name_col, category) triple per row. Dedups across
        the whole ingest via `seen` (first occurrence wins)."""
        n = 0
        batch: list[dict] = []
        for row in self._read_csv(path):
            for id_col, name_col, cat_spec in id_name_cats:
                local = (row.get(id_col) or "").strip()
                if not local:
                    continue
                curie = self._normalize_curie(local)
                if not curie or curie in seen:
                    continue
                seen.add(curie)
                name = (row.get(name_col) or "").strip() or None
                synonyms = {s for c in syn_cols if (s := (row.get(c) or "").strip())}
                category = cat_spec(row) if callable(cat_spec) else cat_spec
                batch.append(self._make_node(curie, name, synonyms, category))
                n += 1
            if len(batch) >= BATCH_SIZE:
                save_to_jsonl(batch, nodes_output, mode="a")
                batch = []
        if batch:
            save_to_jsonl(batch, nodes_output, mode="a")
        return n

    def _normalize_curie(self, local: str) -> str | None:
        """Form the canonical CURIE for a LOINC local id via biomapper2 -- the single source of CURIE/prefix
        logic -- rather than string-building it here, so LOINC prefix/format handling lives in one place and
        stays correct if it ever changes. Returns None (and counts it) if biomapper2 can't validate the id."""
        if local in self._curie_cache:
            return self._curie_cache[local]
        resolved, _, _ = self.normalizer.get_curies(
            {LOINC_PREFIX: local}, stop_on_invalid_id=False, log_warnings=False
        )
        curie = next(iter(resolved), None)
        if not curie:
            self._n_unmapped += 1
            logging.warning(f"{self.source_name}: biomapper2 could not form a CURIE for {LOINC_PREFIX}:{local}; dropping it.")
        self._curie_cache[local] = curie
        return curie

    def _make_node(self, curie: str, name: str | None, synonyms: set[str], category: str) -> dict:
        return self.create_node(
            curie=curie,
            categories=[category],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=name,
            synonyms=(synonyms or None),
        )

    @staticmethod
    def _read_csv(path: Path):
        with open(path, encoding="utf-8") as f:
            yield from csv.DictReader(f)
