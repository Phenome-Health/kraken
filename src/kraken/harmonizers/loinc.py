# loinc.py
import csv
import logging
import re
from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.name_overrides import ORIGINAL_NAME_ATTRIBUTE
from kraken.utils.constants import KNOWLEDGE_ASSERTION, LOINC_INFORES, MANUAL_AGENT
from kraken.utils.kg_io import save_to_jsonl

# LOINC terminology ingest: one node per LOINC identifier, with a display name + name synonyms. Everything
# lives under the single LOINC: prefix (numeric codes, LP parts, LA answers, LL answer lists, LG groups).
# Edges come only from LOINC's own external-terminology mappings (see EXT_CODE_SYSTEMS); LOINC's internal
# hierarchy is not ingested.
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


# --- External code mappings ---
# LOINC maps its Parts to external terminologies (PartRelatedCodeMapping.csv) and its answers to external
# concepts (AnswerList.csv). Both name the target system by URI; these are the ones biomapper2 can form a
# CURIE for. Anything else is counted and logged rather than dropped silently -- RadLex (~1.1k part mappings)
# and ClinVar (~87) are the notable gaps, and both are biomapper2 vocab candidates.
EXT_CODE_SYSTEMS = {
    "http://snomed.info/sct": "SNOMEDCT",
    "https://www.ebi.ac.uk/chebi": "CHEBI",
    "https://www.ncbi.nlm.nih.gov/gene": "NCBIGene",
    "http://www.genenames.org": "HGNC",
    "http://www.nlm.nih.gov/research/umls/rxnorm": "RXCUI",
    "https://www.ncbi.nlm.nih.gov/taxonomy": "NCBITaxon",
    "http://fdasis.nlm.nih.gov": "UNII",
    "http://pubchem.ncbi.nlm.nih.gov": "PUBCHEM.COMPOUND",
}

PART_MAPPING_FILE = "AccessoryFiles/PartFile/PartRelatedCodeMapping.csv"
ANSWER_FILE_PATH = "AccessoryFiles/PanelsAndForms/AnswerList.csv"

# LOINC's `Equivalence` values describe the TARGET concept relative to the LOINC one ("maps in which the
# target concept is more granular than the LOINC concept ... are represented as 'narrower'" -- PartFileReadMe),
# which is the same direction Biolink's match predicates describe their object. So the mapping is direct, with
# LOINC's part as subject.
#
# Deliberately NOT subclass_of: these are cross-terminology mappings, not ontological subsumption. SKOS keeps
# that distinction (`broader` within a scheme vs `broadMatch` across schemes) and so does entity resolution,
# which must not treat a narrow/broad match as an equivalence.
EQUIVALENT = "equivalent"  # becomes an equivalent_id rather than an edge
EQUIVALENCE_PREDICATES = {
    "narrower": "biolink:narrow_match",
    "wider": "biolink:broad_match",
    "relatedto": "biolink:related_to",
}

# An answer code is an information artifact that REFERS to a concept -- "Pompe disease [answer]" is a
# permissible response, not the disorder -- so it mentions its external code rather than equating to it.
# (Those two are also different Biolink branches, so equating them would be an entity-resolution error.)
ANSWER_REFERENCE_PREDICATE = "biolink:mentions"


@dataclass(frozen=True)
class NodeSpec:
    """One kind of node to mint from a row. `category` and `description` are each either a fixed value, a
    column name (description only), or a callable(row) for per-row logic."""

    id_col: str
    name_col: str
    category: str | Callable[[dict], str]
    description: str | Callable[[dict], str | None] | None = None
    # Columns a description CALLABLE reads. Needed so they aren't also retained verbatim in attributes; when
    # `description` is a plain column name that's inferred instead.
    description_cols: tuple[str, ...] = ()


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
        (
            NodeSpec(
                "LOINC_NUM",
                "LONG_COMMON_NAME",
                numeric_category,
                numeric_description,
                description_cols=(DEFINITION_COLUMN, SURVEY_QUESTION_COLUMN),
            ),
        ),
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
        ANSWER_FILE_PATH,
        (
            NodeSpec("AnswerStringId", "DisplayText", INFO_CATEGORY, "Description"),
            NodeSpec("AnswerListId", "AnswerListName", INFO_CATEGORY),
        ),
    ),
]

BATCH_SIZE = 50_000


class LoincHarmonizer(BaseHarmonizer):
    """LOINC ingest. Mints one node per LOINC identifier (numeric lab codes, LP parts, LA answers, LL answer
    lists, LG groups) with a display name and any short/display-name synonyms, all under the LOINC: prefix.

    Takes LOINC's own mappings to external terminologies: `equivalent` ones become equivalent_ids (giving
    otherwise-orphan Parts a route into the graph), and the narrower/wider/relatedto ones become edges, as do
    the concepts an answer code refers to. LOINC's internal hierarchy is NOT ingested.
    """

    source_infores = LOINC_INFORES

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self._curie_cache: dict[str, str | None] = {}  # LOINC local id -> canonical curie (or None)
        self._n_unmapped = 0
        self._n_qualified = 0
        self._edges: list[dict] = []
        self._part_equivalences: dict[str, list[str]] = {}
        self._unmapped_code_systems: dict[str, int] = defaultdict(int)
        self._invalid_ext_codes: dict[str, int] = defaultdict(int)
        self._unknown_equivalences: dict[str, int] = defaultdict(int)
        self._external_nodes: dict[str, str] = {}  # edge target CURIE -> its display name

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

        self._n_unmapped = 0
        self._n_qualified = 0
        seen: set[str] = set()
        counts: dict[str, int] = {}
        # External mappings first: the equivalent ones become equivalent_ids on the part nodes minted below,
        # and the rest (plus the answers' references) become this source's only edges.
        self._part_equivalences = self._load_part_mappings(base)
        self._collect_answer_reference_edges(base)

        for loinc_file in LOINC_FILES:
            counts[loinc_file.path] = self._ingest(base, loinc_file, nodes_output, seen)
        counts["external stubs"] = self._write_external_nodes(nodes_output, seen)
        save_to_jsonl(self._edges, edges_output, mode="w")

        summary = ", ".join(f"{Path(p).name}={n}" for p, n in counts.items())
        logging.info(
            f"{self.source_name} harmonization complete: {len(seen)} LOINC nodes, {len(self._edges)} edges "
            f"({summary})"
            + (f"; {self._n_qualified} part names qualified by part type" if self._n_qualified else "")
            + (f"; {self._n_unmapped} ids dropped (biomapper2 could not form a CURIE)" if self._n_unmapped else "")
        )
        if self._unmapped_code_systems:
            logging.warning(
                f"Dropped external mappings from {len(self._unmapped_code_systems)} code system(s) with no "
                f"biomapper2 vocab: {dict(self._unmapped_code_systems)}. Register the worthwhile ones there."
            )
        if self._invalid_ext_codes:
            logging.warning(f"External codes that failed biomapper2 validation: {dict(self._invalid_ext_codes)}")
        if self._unknown_equivalences:
            logging.warning(
                f"Unrecognized Equivalence values (mapping neither kept nor edged): "
                f"{dict(self._unknown_equivalences)}"
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
                # Retain every column we didn't map to a property, matching what BaseHarmonizer does for the
                # generic sources. LOINC's six axes (COMPONENT, PROPERTY, TIME_ASPCT, SYSTEM, SCALE_TYP,
                # METHOD_TYP) are what a code actually IS, and they'd otherwise be dropped.
                attributes |= self._unmapped_columns(row, loinc_file, spec)

                batch.append(
                    self._make_node(
                        curie,
                        name,
                        synonyms,
                        category,
                        description,
                        attributes,
                        self._part_equivalences.get(curie, ()),
                    )
                )
                n += 1
            if len(batch) >= BATCH_SIZE:
                save_to_jsonl(batch, nodes_output, mode="a")
                batch = []
        if batch:
            save_to_jsonl(batch, nodes_output, mode="a")
        return n

    # ------------------------------------------------------------------ external code mappings

    def _ext_curie(self, row: dict) -> str | None:
        """CURIE for a row's external code, or None if we can't form one. LOINC names the target system by
        URI and sometimes writes the code already prefixed ("CHEBI:33216"), so the prefix is stripped before
        biomapper2 -- which owns what a valid id looks like -- is asked to build the CURIE."""
        system = (row.get("ExtCodeSystem") or "").strip()
        code = (row.get("ExtCodeId") or "").strip()
        if not code:
            return None
        vocab = EXT_CODE_SYSTEMS.get(system)
        if not vocab:
            self._unmapped_code_systems[system or "(blank)"] += 1
            return None
        local_id = code.split(":", 1)[1] if code.upper().startswith(f"{vocab.upper()}:") else code
        resolved, invalid, _ = self.normalizer.get_curies(
            {vocab: local_id}, stop_on_invalid_id=False, log_warnings=False, fuzzy_match_vocab=False
        )
        if invalid:
            self._invalid_ext_codes[vocab] += 1
        return next(iter(resolved), None)

    def _load_part_mappings(self, base: Path) -> dict[str, list[str]]:
        """Read PartRelatedCodeMapping, returning part CURIE -> equivalent external CURIEs and stashing the
        non-equivalent mappings as edges. LOINC states the specificity of each mapping, so `equivalent` can be
        trusted as an equivalence while narrower/wider/relatedto must stay edges."""
        path = base / PART_MAPPING_FILE
        if not path.is_file():
            logging.warning(f"{self.source_name}: no {PART_MAPPING_FILE}; parts will have no external mappings")
            return {}

        equivalences: dict[str, list[str]] = defaultdict(list)
        for row in self._read_csv(path):
            part_curie = self._normalize_curie((row.get("PartNumber") or "").strip())
            ext_curie = self._ext_curie(row)
            if not part_curie or not ext_curie:
                continue
            equivalence = (row.get("Equivalence") or "").strip().lower()
            if equivalence == EQUIVALENT:
                equivalences[part_curie].append(ext_curie)
            elif predicate := EQUIVALENCE_PREDICATES.get(equivalence):
                self._edges.append(self._make_edge(part_curie, predicate, ext_curie, row))
            else:
                self._unknown_equivalences[equivalence or "(blank)"] += 1
        logging.info(
            f"Loaded {sum(len(v) for v in equivalences.values())} equivalent external codes for "
            f"{len(equivalences)} LOINC parts, plus {len(self._edges)} non-equivalent mapping edges"
        )
        self._warn_about_shared_equivalences(equivalences)
        return equivalences

    def _warn_about_shared_equivalences(self, equivalences: dict[str, list[str]]):
        """Report parts that share an equivalent external code, since those will merge transitively.

        LOINC's mappings are only as fine-grained as the terminology they point at: Phenytoin.free, .bound and
        .total all map to one PubChem code, as do Triiodothyronine and its reverse isomer. We take LOINC at its
        word anyway -- it's the source's own assertion, and no local rule reliably separates a good mapping
        from a coarse one without knowing the target's granularity -- but entity resolution should weight these
        below corroborated equivalences, and this is the number that says how much is at stake."""
        claimants: dict[str, list[str]] = defaultdict(list)
        for part_curie, ext_curies in equivalences.items():
            for ext_curie in ext_curies:
                claimants[ext_curie].append(part_curie)
        shared = {code: parts for code, parts in claimants.items() if len(parts) > 1}
        if not shared:
            return
        affected = {part for parts in shared.values() for part in parts}
        logging.warning(
            f"{len(shared)} external codes are claimed as equivalent by more than one LOINC part, pulling "
            f"{len(affected)} parts into transitive merges (largest group {max(len(p) for p in shared.values())}). "
            f"LOINC's mappings are as coarse as the target terminology, so these are a known over-merge risk "
            f"and warrant a lower weight in entity resolution."
        )

    def _collect_answer_reference_edges(self, base: Path):
        """An answer's external code says which concept it refers to, which is a `mentions` edge rather than
        an equivalence (see ANSWER_REFERENCE_PREDICATE)."""
        path = base / ANSWER_FILE_PATH
        if not path.is_file():
            return
        seen_pairs: set[tuple[str, str]] = set()  # an answer repeats across the lists it belongs to
        for row in self._read_csv(path):
            answer_curie = self._normalize_curie((row.get("AnswerStringId") or "").strip())
            ext_curie = self._ext_curie(row)
            if not answer_curie or not ext_curie or (answer_curie, ext_curie) in seen_pairs:
                continue
            seen_pairs.add((answer_curie, ext_curie))
            self._edges.append(self._make_edge(answer_curie, ANSWER_REFERENCE_PREDICATE, ext_curie, row))

    def _write_external_nodes(self, nodes_output: Path, seen: set[str]) -> int:
        """Mint a stub node for each external concept our edges point at, since every edge endpoint must exist
        as a node in the same artifact.

        Typed NamedThing: we don't know what an arbitrary SNOMED or ChEBI concept is, and NamedThing is the
        honest placeholder -- it's also a wildcard for taxon/branch merge guards, so the stub takes on the real
        type when the authoritative node for that CURIE merges in. LOINC does give us a display name for each,
        so these aren't nameless."""
        stubs = [
            self._make_node(curie, display_name or None, set(), NAMED_THING)
            for curie, display_name in sorted(self._external_nodes.items())
            if curie not in seen
        ]
        if stubs:
            save_to_jsonl(stubs, nodes_output, mode="a")
            seen.update(node["id"] for node in stubs)
        return len(stubs)

    def _make_edge(self, subject: str, predicate: str, obj: str, row: dict) -> dict:
        """LOINC curates these mappings by hand, hence knowledge_assertion / manual_agent.

        Also records the edge's target so a stub node can be minted for it: every edge endpoint has to exist
        as a node in the same artifact, and these point at external terminologies KRAKEN may not carry yet."""
        self._external_nodes.setdefault(obj, (row.get("ExtCodeDisplayName") or "").strip())
        return self.create_edge(
            subject_id=subject,
            object_id=obj,
            predicate=predicate,
            primary_ks=self.source_infores,
            knowledge_level=KNOWLEDGE_ASSERTION,
            agent_type=MANUAL_AGENT,
        )

    @staticmethod
    def _unmapped_columns(row: dict, loinc_file: LoincFile, spec: NodeSpec) -> dict:
        """Every non-empty column that didn't become a property, kept verbatim. On the answer file each row
        mints two nodes, so each also carries the other's columns -- which is useful (an answer keeps the id
        of the list it belongs to) rather than noise."""
        consumed = {spec.id_col, spec.name_col, *loinc_file.synonym_cols, *spec.description_cols}
        if isinstance(spec.description, str):
            consumed.add(spec.description)
        return {k: v.strip() for k, v in row.items() if k not in consumed and (v or "").strip()}

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
        external_equivalent_ids: tuple[str, ...] = (),
    ) -> dict:
        return self.create_node(
            curie=curie,
            categories=[category],
            provided_by=self.source_infores,
            equivalent_ids=[curie, *external_equivalent_ids],
            name=name,
            synonyms=(synonyms or None),
            description=description,
            attributes=(attributes or None),
        )

    @staticmethod
    def _read_csv(path: Path):
        with open(path, encoding="utf-8") as f:
            yield from csv.DictReader(f)
