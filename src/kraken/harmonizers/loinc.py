# loinc.py
import csv
import logging
import re
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.name_overrides import ORIGINAL_NAME_ATTRIBUTE
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


# --- LP part types ---
# A LOINC Part is one axis-value of a LOINC term, and its PartTypeName says which axis. The same PartName can
# appear under several types -- "Creatinine" is both a COMPONENT (LP14355-9, the analyte being measured) and a
# DIVISORS (LP32035-5, the denominator of a ratio) -- so a bare name is ambiguous, and 1,711 of LOINC's 72k
# part names are used by more than one type. The type is appended in brackets ("Creatinine [divisor]") and the
# raw value kept in attributes.
#
# These are NOT merged with one another: LOINC mints them as distinct identifiers and asserts no equivalence
# between them, so combining them would be our invention. Real equivalences still flow in from sources that do
# assert them -- the COMPONENT part above already merges into the CHEBI creatinine node.
PART_TYPE_COLUMN = "PartTypeName"
# COMPONENT is LOINC's unmarked default and 85% of all parts; a component part simply IS the analyte, so
# qualifying those would put noise on 63k nodes for no gain.
UNQUALIFIED_PART_TYPE = "COMPONENT"
# LOINC's own values vary in case and number; displayed lowercase and singular. Only the plural needs stating
# (naive de-pluralizing would turn CLASS into "clas").
PART_TYPE_DISPLAY_OVERRIDES = {"DIVISORS": "divisor"}
# Splits camelCase values like "SubjectMatterDomain" -> "Subject Matter Domain"
_CAMEL_CASE_BOUNDARY = re.compile(r"(?<=[a-z])(?=[A-Z])")


def part_type_display(part_type: str) -> str | None:
    """The bracketed term for a part type, or None when the part needs no qualifier (COMPONENT, or blank).

    Values are lowercased, and dotted ones reduced to their FIRST and LAST segments: the first names the
    family the axis belongs to, the last is the axis itself, and any middle segments are elaboration that
    reads as clutter. "Document.Kind" -> "document kind" (the family is what makes "kind" mean anything),
    "Rad.Guidance for.Action" -> "rad action", "Rad.Modality.Modality Subtype" -> "rad modality subtype"
    (which also avoids the "modality modality" a naive join would produce)."""
    part_type = (part_type or "").strip()
    if not part_type or part_type.upper() == UNQUALIFIED_PART_TYPE:
        return None
    if part_type in PART_TYPE_DISPLAY_OVERRIDES:
        return PART_TYPE_DISPLAY_OVERRIDES[part_type]

    segments = part_type.split(".")
    term = segments[0] if len(segments) == 1 else f"{segments[0]} {segments[-1]}"
    return _CAMEL_CASE_BOUNDARY.sub(" ", term).lower()


def part_name_qualifier(row: dict) -> str | None:
    """The bracketed term to append to a part's name, or None to leave it alone. Skipped when the name already
    says what it is -- LOINC's GENE parts are all named "<symbol> gene", and many SYSTEM parts contain
    "system" -- so this is idempotent and safe to re-run."""
    display = part_type_display(row.get(PART_TYPE_COLUMN, ""))
    if display and display in (row.get("PartName") or "").lower():
        return None
    return display


# --- Descriptions ---
# LOINC writes a definition for only a minority of its codes, and none at all for Parts (Part.csv has no such
# column). Where a numeric code has no definition but is a survey item, the question as actually asked is the
# best description available -- and worth taking even though it often echoes the name, since 60% of the time
# it differs, sometimes substantively ("Race during assessment period [CMS Assessment]" is asked as "What is
# your race?", "Breathing" as "Breathing independent of vocalization").
DEFINITION_COLUMN = "DefinitionDescription"
SURVEY_QUESTION_COLUMN = "SURVEY_QUEST_TEXT"
SURVEY_DESCRIPTION_PREFIX = "Survey question: "


def numeric_description(row: dict) -> str | None:
    """Description for a numeric LOINC code: its definition, else its survey question text (labelled, so the
    description says what kind of thing it is rather than just restating the name)."""
    if definition := (row.get(DEFINITION_COLUMN) or "").strip():
        return definition
    if question := (row.get(SURVEY_QUESTION_COLUMN) or "").strip():
        return f"{SURVEY_DESCRIPTION_PREFIX}{question}"
    return None


@dataclass(frozen=True)
class NodeSpec:
    """One kind of node to mint from a row. `category` and `description` are each either a fixed value, a
    column name (description only), or a callable(row) for per-row logic."""

    id_col: str
    name_col: str
    category: str | Callable[[dict], str]
    description: str | Callable[[dict], str | None] | None = None


@dataclass(frozen=True)
class LoincFile:
    """A LOINC CSV and how to read nodes out of it. Most yield one node per row; AnswerList.csv yields two
    (the answer and the answer list it belongs to)."""

    path: str
    nodes: tuple[NodeSpec, ...]
    synonym_cols: tuple[str, ...] = ()
    name_qualifier: Callable[[dict], str | None] | None = None


LOINC_FILES = [
    LoincFile(
        "LoincTable/Loinc.csv",
        (NodeSpec("LOINC_NUM", "LONG_COMMON_NAME", numeric_category, numeric_description),),
        synonym_cols=("SHORTNAME", "DisplayName", "CONSUMER_NAME"),
    ),
    LoincFile(
        "AccessoryFiles/PartFile/Part.csv",
        (NodeSpec("PartNumber", "PartName", NAMED_THING),),
        synonym_cols=("PartDisplayName",),
        name_qualifier=part_name_qualifier,
    ),
    LoincFile("AccessoryFiles/GroupFile/Group.csv", (NodeSpec("GroupId", "Group", INFO_CATEGORY),)),
    LoincFile("AccessoryFiles/GroupFile/ParentGroup.csv", (NodeSpec("ParentGroupId", "ParentGroup", INFO_CATEGORY),)),
    # Each row carries both an answer (LA...) and the answer list it belongs to (LL...). The Description
    # column describes the answer, so only that spec takes it.
    LoincFile(
        "AccessoryFiles/PanelsAndForms/AnswerList.csv",
        (
            NodeSpec("AnswerStringId", "DisplayText", INFO_CATEGORY, "Description"),
            NodeSpec("AnswerListId", "AnswerListName", INFO_CATEGORY),
        ),
    ),
]

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
        self._n_qualified = 0

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
        self._n_qualified = 0
        seen: set[str] = set()
        counts: dict[str, int] = {}
        for loinc_file in LOINC_FILES:
            counts[loinc_file.path] = self._ingest(base, loinc_file, nodes_output, seen)

        summary = ", ".join(f"{Path(p).name}={n}" for p, n in counts.items())
        logging.info(
            f"{self.source_name} harmonization complete: {len(seen)} LOINC nodes ({summary})"
            + (f"; {self._n_qualified} part names qualified by part type" if self._n_qualified else "")
            + (f"; {self._n_unmapped} ids dropped (biomapper2 could not form a CURIE)" if self._n_unmapped else "")
        )

    # ------------------------------------------------------------------ core

    def _ingest(self, base: Path, loinc_file: LoincFile, nodes_output: Path, seen: set[str]) -> int:
        """Stream a LOINC CSV, minting a node for each of the file's NodeSpecs per row. Dedups across the whole
        ingest via `seen` (first occurrence wins)."""
        n = 0
        batch: list[dict] = []
        for row in self._read_csv(base / loinc_file.path):
            for spec in loinc_file.nodes:
                local = (row.get(spec.id_col) or "").strip()
                if not local:
                    continue
                curie = self._normalize_curie(local)
                if not curie or curie in seen:
                    continue
                seen.add(curie)
                name = (row.get(spec.name_col) or "").strip() or None
                synonyms = {s for c in loinc_file.synonym_cols if (s := (row.get(c) or "").strip())}
                category = spec.category(row) if callable(spec.category) else spec.category
                description = self._description(row, spec)

                # Say which axis a part belongs to, so "Creatinine" the divisor is distinguishable from
                # "Creatinine" the analyte. The raw LOINC value is kept for anyone who needs it back.
                attributes = {}
                qualifier = loinc_file.name_qualifier(row) if loinc_file.name_qualifier else None
                if qualifier and name:
                    # Same attribute key the cross-source name overrides use (harmonizers/helpers/
                    # name_overrides.py), so "what was this called before we qualified it?" has one answer
                    # however the node was renamed. Synonyms aren't stored separately -- they carry the same
                    # bracket, so stripping it recovers them.
                    attributes[ORIGINAL_NAME_ATTRIBUTE] = name
                    name = self._qualify(name, qualifier)
                    self._n_qualified += 1
                    # Qualify the synonyms too. LOINC repeats PartName verbatim as PartDisplayName for 73% of
                    # parts, so leaving them bare would put an unqualified "Creatinine" back on the divisor
                    # node -- exactly the collision with the analyte that qualifying is meant to remove.
                    synonyms = {self._qualify(synonym, qualifier) for synonym in synonyms}
                if raw_part_type := (row.get(PART_TYPE_COLUMN) or "").strip():
                    attributes[PART_TYPE_COLUMN] = raw_part_type

                batch.append(self._make_node(curie, name, synonyms, category, description, attributes))
                n += 1
            if len(batch) >= BATCH_SIZE:
                save_to_jsonl(batch, nodes_output, mode="a")
                batch = []
        if batch:
            save_to_jsonl(batch, nodes_output, mode="a")
        return n

    @staticmethod
    def _description(row: dict, spec: NodeSpec) -> str | None:
        """A spec's description: nothing, a named column, or a callable that picks per row."""
        if spec.description is None:
            return None
        if callable(spec.description):
            return spec.description(row)
        return (row.get(spec.description) or "").strip() or None

    @staticmethod
    def _qualify(text: str, qualifier: str) -> str:
        """Append the bracketed part type, unless the text already says it (so this stays idempotent)."""
        return text if qualifier in text.lower() else f"{text} [{qualifier}]"

    def _normalize_curie(self, local: str) -> str | None:
        """Form the canonical CURIE for a LOINC local id via biomapper2 -- the single source of CURIE/prefix
        logic -- rather than string-building it here, so LOINC prefix/format handling lives in one place and
        stays correct if it ever changes. Returns None (and counts it) if biomapper2 can't validate the id."""
        if local in self._curie_cache:
            return self._curie_cache[local]
        resolved, _, _ = self.normalizer.get_curies({LOINC_PREFIX: local}, stop_on_invalid_id=False, log_warnings=False)
        curie = next(iter(resolved), None)
        if not curie:
            self._n_unmapped += 1
            logging.warning(
                f"{self.source_name}: biomapper2 could not form a CURIE for {LOINC_PREFIX}:{local}; dropping it."
            )
        self._curie_cache[local] = curie
        return curie

    def _make_node(
        self,
        curie: str,
        name: str | None,
        synonyms: set[str],
        category: str,
        description: str | None = None,
        attributes: dict | None = None,
    ) -> dict:
        return self.create_node(
            curie=curie,
            categories=[category],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=name,
            synonyms=(synonyms or None),
            description=description,
            attributes=(attributes or None),
        )

    @staticmethod
    def _read_csv(path: Path):
        with open(path, encoding="utf-8") as f:
            yield from csv.DictReader(f)
