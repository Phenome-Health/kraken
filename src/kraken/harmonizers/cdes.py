# cdes.py
import json
import logging
from pathlib import Path
from typing import Any

from kraken.biolink_client import BiolinkClient
from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.cde_concept_blocklist import CONCEPT_BLOCKLIST
from kraken.utils.constants import KNOWLEDGE_ASSERTION, MANUAL_AGENT, NIH_CDE_SOURCE_ID
from kraken.utils.general import clean_text
from kraken.utils.kg_io import save_to_jsonl

# NIH CDE source vocab name (as it appears in the export) -> biomapper2 vocab name.
# Note: we intentionally exclude SNOMEDCT, due to licensing restrictions.
VOCAB_MAP = {
    "nci thesaurus": "ncit",
    "ncit": "ncit",
    "nci": "ncit",
    "umls": "umls",
    "loinc": "loinc",
    "lnc": "loinc",
}

# caDSR-origin concepts are registry composites (not ontology terms) -> attributes only, never edges
CADSR_ORIGINS = {"nci cadsr", "cadsr"}

# Priority for choosing the CDE node's primary name; the rest become synonyms. We favor the
# formal element *name* over survey *question wording* -- e.g. "Address City Name" over the
# Preferred Question Text "City" -- since the question phrasing is better kept as a synonym.
# UNTAGGED_NAME is a sentinel for the primary (untagged) designation, which is usually the
# formal name (e.g. "Address City Name") when no Long/Full Name tag is present.
UNTAGGED_NAME = "<untagged>"
NAME_TAG_PRIORITY = [
    "Long Common Name",  # LOINC-style; usually clean and readable
    "Full Name",
    UNTAGGED_NAME,  # primary designation, e.g. "Address City Name"
    "Preferred Question Text",
    "Question Text",
    "Suggested Question Text",
    "Long Name",  # caDSR-style; often machine-concatenated, so below question wording
    "Short Name",
    "Shortname",
    "Display Name",
]

# Words too generic to stand alone as a CDE name (structural/metadata terms & response glue). When
# the highest-priority designation is exactly one of these, we fall back to the next, more specific
# designation -- e.g. untagged "Etiology" -> "Etiology liver abnormality unknown indicate code".
# Deliberately excludes substantive standalone concepts (Age, Race, Sex, disease/drug/lab names).
GENERIC_NAMES: frozenset[str] = frozenset({
    "etiology", "diagnosis", "type", "status", "code", "score", "result", "indicator", "category",
    "criterion", "value", "level", "grade", "stage", "class", "group", "method", "measure",
    "outcome", "reason", "response", "rating", "unit", "amount", "count", "number", "frequency",
    "duration", "date", "time", "name", "description", "comment", "comments", "other", "specify",
    "describe", "present", "protocol", "center", "site", "source", "location", "organism",
})

CDE_PREFIX = "CDE"
CDE_CATEGORY = "biolink:CommonDataElement"
STUB_CATEGORY = "biolink:NamedThing"
ASSESSES = "biolink:assesses"  # CDE -> a concept it measures
RELATED_TO = "biolink:related_to"  # CDE -> a permissible-answer concept


class CDEHarmonizer(BaseHarmonizer):
    """Harmonizer for the NIH CDE Repository export (a single nested JSON array of CDEs).

    Each CDE becomes a ``biolink:CommonDataElement`` node. LOINC ids go into ``equivalent_ids``
    (so the CDE merges into / names kraken's LOINC nodes). Ontology concept mappings
    (objectClass / property / dataElementConcept) become ``assesses`` edges, and coded permissible
    values become ``related_to`` edges -- both filtered by a curated concept blocklist. Compound
    codes (e.g. "Myocarditis/Pericarditis") are split into their constituent concepts by the
    normalizer. A minimal, named stub node is minted for every edge endpoint so nothing orphans;
    stubs carry only ``biolink:NamedThing`` and get properly typed when richer sources merge.
    """

    source_name = "cdes"
    # Not a registered infores; used verbatim as the source id per project decision.
    source_infores = NIH_CDE_SOURCE_ID

    def __init__(self, biolink_client: BiolinkClient):
        super().__init__(biolink_client)
        self._stub_nodes: dict[str, dict] = {}  # curie -> minimal node (deduped across CDEs)
        self._curie_cache: dict[tuple, list[str]] = {}  # (vocab, code) -> resolved curie(s)

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
            raise ValueError(f"{self.source_name} requires input_file")
        logging.info(f"Harmonizing {self.source_name}: {input_file} -> {nodes_output}, {edges_output}")

        with open(input_file) as f:
            records = json.load(f)

        cde_nodes, edges = [], []
        for record in records:
            node, record_edges = self._harmonize_cde(record)
            if node:
                cde_nodes.append(node)
                edges.extend(record_edges)

        all_nodes = cde_nodes + list(self._stub_nodes.values())
        save_to_jsonl(all_nodes, nodes_output, mode="w")
        save_to_jsonl(edges, edges_output, mode="w")
        logging.info(
            f"{self.source_name} harmonization complete: {len(cde_nodes)} CDE nodes + "
            f"{len(self._stub_nodes)} stub nodes, {len(edges)} edges"
        )

    # ------------------------------------------------------------------ per-record

    def _harmonize_cde(self, record: dict[str, Any]) -> tuple[dict | None, list[dict]]:
        tiny_id = record.get("tinyId")
        if not tiny_id:
            return None, []
        cde_curie = f"{CDE_PREFIX}:{tiny_id}"

        name, synonyms = self._select_name_and_synonyms(record.get("designations", []))
        description = self._select_description(record.get("definitions", []))

        # LOINC ids -> equivalent_ids (merge into + enrich kraken's LOINC nodes)
        equivalent_ids = [cde_curie]
        for idn in record.get("ids", []):
            if (idn.get("source") or "").strip().lower() == "loinc":
                equivalent_ids.extend(self._construct_curies("loinc", idn.get("id")))

        cde_node = self.create_node(
            curie=cde_curie,
            categories=[CDE_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=equivalent_ids,
            name=name,
            synonyms=synonyms,
            description=description,
            attributes=self._collect_cde_attributes(record),
        )

        edges: list[dict] = []
        seen_objects: set[str] = set()

        # assesses edges: ontology-origin concepts only (caDSR composites stay in attributes)
        for field in ("objectClass", "property", "dataElementConcept"):
            for concept in record.get(field, {}).get("concepts", []):
                for curie in self._construct_curies(concept.get("origin"), concept.get("originId")):
                    if self._edge_worthy(curie, seen_objects):
                        seen_objects.add(curie)
                        self._add_stub(curie, concept.get("name"))
                        edges.append(self._make_edge(cde_curie, curie, ASSESSES))

        # related_to edges: coded permissible values. Compound codes (e.g. "Myocarditis/Pericarditis")
        # split into multiple concepts; only single-concept answers carry the answer label onto the stub.
        for pv in record.get("valueDomain", {}).get("permissibleValues", []):
            answer_curies = self._answer_curies(pv)
            label = (pv.get("valueMeaningName") or pv.get("permissibleValue")) if len(answer_curies) == 1 else None
            for curie in answer_curies:
                if self._edge_worthy(curie, seen_objects):
                    seen_objects.add(curie)
                    self._add_stub(curie, label)
                    edges.append(self._make_edge(cde_curie, curie, RELATED_TO))

        return cde_node, edges

    @staticmethod
    def _edge_worthy(curie: str | None, already_seen: set[str]) -> bool:
        return bool(curie) and curie not in CONCEPT_BLOCKLIST and curie not in already_seen

    # ------------------------------------------------------------------ curie construction

    def _construct_curies(self, system: str | None, code: str | None) -> list[str]:
        """Map a source vocab name + local code to validated CURIE(s).

        Compound codes (e.g. "C34831:C34915" for a slash-answer like "Myocarditis/Pericarditis")
        are split on ':' by the normalizer and each code resolved independently; unparseable
        fragments and label-as-code non-codes ("N/A or not reported") drop out via validation.
        Returns [] for unrecognized vocabs or when nothing resolves. All id handling lives in
        biomapper2 -- we just pass the ':' delimiter.
        """
        if not system or not code:
            return []
        code = str(code).strip()
        if not code:
            return []
        vocab = VOCAB_MAP.get(system.strip().lower())
        if not vocab:
            return []

        cache_key = (vocab, code)
        if cache_key not in self._curie_cache:
            curies, _, _ = self.normalizer.get_curies(
                {vocab: code},
                stop_on_invalid_id=False,
                log_warnings=False,
                fuzzy_match_vocab=False,
                array_delimiters=[":"],
            )
            self._curie_cache[cache_key] = list(curies)
        return self._curie_cache[cache_key]

    def _answer_curies(self, pv: dict[str, Any]) -> list[str]:
        # Prefer the concept slot (usually clean NCIt), fall back to the code slot (usually LOINC answers)
        return self._construct_curies(pv.get("conceptSource"), pv.get("conceptId")) or self._construct_curies(
            pv.get("codeSystemName"), pv.get("valueMeaningCode")
        )

    # ------------------------------------------------------------------ node/edge builders

    def _add_stub(self, curie: str, name: str | None):
        if curie not in self._stub_nodes:
            self._stub_nodes[curie] = self.create_node(
                curie=curie,
                categories=[STUB_CATEGORY],
                provided_by=self.source_infores,
                equivalent_ids=[curie],
                name=name,
            )

    def _make_edge(self, subject: str, obj: str, predicate: str) -> dict:
        return self.create_edge(
            subject_id=subject,
            object_id=obj,
            predicate=predicate,
            primary_ks=self.source_infores,
            knowledge_level=KNOWLEDGE_ASSERTION,
            agent_type=MANUAL_AGENT,
        )

    # ------------------------------------------------------------------ field selection

    def _select_name_and_synonyms(self, designations: list[dict]) -> tuple[str | None, list[str]]:
        by_tag: dict[str, str] = {}
        all_texts: list[str] = []
        for d in designations:
            text = (d.get("designation") or "").strip()  # create_node's clean_text collapses inner whitespace
            if not text:
                continue
            all_texts.append(text)
            if d.get("tags"):
                for tag in d["tags"]:
                    by_tag.setdefault(tag, text)
            else:
                by_tag.setdefault(UNTAGGED_NAME, text)
        # Ordered name candidates by priority, with the first designation as a last resort.
        candidates = [by_tag[t] for t in NAME_TAG_PRIORITY if t in by_tag]
        if all_texts:
            candidates.append(all_texts[0])
        name = candidates[0] if candidates else None
        # If the top pick is a too-generic standalone word, fall back to the first more-specific
        # (longer, non-generic) designation -- but keep the generic word if the only alternatives are
        # shorter abbreviations (e.g. "Diagnosis" is a better name than "Dx").
        if name and name.strip().lower() in GENERIC_NAMES:
            better = next(
                (c for c in candidates if len(c) > len(name) and c.strip().lower() not in GENERIC_NAMES),
                None,
            )
            name = better or name
        synonyms = [t for t in all_texts if t != name]
        return name, synonyms

    @staticmethod
    def _select_description(definitions: list[dict]) -> str | None:
        # clean_text here so the longest-pick compares collapsed text (create_node cleans it again, idempotently)
        raw = [d.get("definition") for d in definitions]
        texts = [clean_text(t) for t in raw if isinstance(t, str) and t.strip()]
        return max(texts, key=len) if texts else None  # longest avoids truncated fragments

    def _collect_cde_attributes(self, record: dict[str, Any]) -> dict[str, Any]:
        vd = record.get("valueDomain", {})
        reg = record.get("registrationState", {})
        attrs: dict[str, Any] = {
            "steward": record.get("steward"),
            "datatype": vd.get("datatype"),
            "registration_status": reg.get("registrationStatus"),
            "administrative_status": reg.get("administrativeStatus"),
            "nih_endorsed": record.get("nihEndorsed"),
            "value_domain_name": vd.get("name"),
            "views": record.get("views"),
        }

        source_ids = [
            {k: v for k, v in {"source": i.get("source"), "id": i.get("id"), "version": i.get("version")}.items() if v}
            for i in record.get("ids", [])
        ]
        if source_ids:
            attrs["source_ids"] = source_ids

        cadsr = [
            {"name": c.get("name"), "id": c.get("originId")}
            for c in record.get("dataElementConcept", {}).get("concepts", [])
            if (c.get("origin") or "").strip().lower() in CADSR_ORIGINS
        ]
        if cadsr:
            attrs["cadsr_data_element_concepts"] = cadsr

        pvs = []
        for pv in vd.get("permissibleValues", []):
            entry = {
                "value": pv.get("permissibleValue"),
                "label": pv.get("valueMeaningName"),
                "code": pv.get("valueMeaningCode"),
                "code_system": pv.get("codeSystemName"),
                "concept_id": pv.get("conceptId"),
                "concept_source": pv.get("conceptSource"),
            }
            entry = {k: v for k, v in entry.items() if v}
            if entry:
                pvs.append(entry)
        if pvs:
            attrs["permissible_values"] = pvs

        # Drop empties (booleans/0 are meaningful and kept)
        return {k: v for k, v in attrs.items() if v is not None and v != "" and v != [] and v != {}}
