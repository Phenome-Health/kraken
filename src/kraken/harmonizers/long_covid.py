# long_covid.py
import csv
import logging
import re
from pathlib import Path
from typing import Any

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import (
    DATA_ANALYSIS_PIPELINE,
    ROOT_CATEGORY,
    STATISTICAL_ASSOCIATION,
)
from kraken.utils.kg_io import save_to_jsonl

# --- Biolink types (each a single swappable constant; approved by Qi Wei, ISB) ---
HAS_PHENOTYPE = "biolink:has_phenotype"  # disease -> symptom (SPOKE DpS)
POSITIVELY_CORRELATED = "biolink:positively_correlated_with"  # disease -> protein (SPOKE INCREASEDIN_PiD)
NEGATIVELY_CORRELATED = "biolink:negatively_correlated_with"  # disease -> protein (SPOKE DECREASEDIN_PdD)
APPLIED_TO_TREAT = "biolink:applied_to_treat"  # compound -> disease
DISEASE_CATEGORY = "biolink:Disease"
PROTEIN_CATEGORY = "biolink:Protein"
COMPOUND_CATEGORY = "biolink:ChemicalEntity"
SYMPTOM_CATEGORY = "biolink:DiseaseOrPhenotypicFeature"  # the symptom tables are mostly phenotypes, some diseases
# Symptom-table MeSH terms that are neither a disease nor a phenotype -- an organ, body fluids, normal functions -- so
# they get the wildcard root category rather than one that would clash with their real type in entity resolution.
# (Presumably the table meant an abnormality of each, e.g. hearing loss; flagged to ISB.) Keyed by symptom_id.
NON_PHENOTYPE_SYMPTOM_IDS = {
    "D001743",  # Urinary Bladder
    "D013183",  # Sputum
    "D013542",  # Sweat
    "D006309",  # Hearing
    "D014785",  # Vision, Ocular
    "D013894",  # Thirst
}

# The vocabulary of each unprefixed id column. The prefixed curie is only a hand-off to create_node/create_edge,
# which run it through biomapper2 for the canonical spelling.
SYMPTOM_VOCAB = "MESH"
PROTEIN_VOCAB = "UniProtKB"
COMPOUND_VOCAB = "CHEMBL.COMPOUND"  # a bare "CHEMBL" vocab resolves to CHEMBL.MECHANISM

# Edge types as SPOKE (and so the protein table's edge_type column) spells them
EDGE_TYPE_TO_PREDICATE = {
    "INCREASEDIN_PiD": POSITIVELY_CORRELATED,
    "DECREASEDIN_PdD": NEGATIVELY_CORRELATED,
}

# Per-table publications (from Qi Wei); the protein table carries its own per-row pmid_list instead
INCOV_PMIDS = ["PMID:35216672"]
INSIGHT_ONEFLORIDA_PMIDS = ["PMID:37029117", "PMID:36785842"]
RECOVER_SYMPTOM_PMIDS = ["PMID:37278994"]

SYMPTOM_FILES = {
    "INCOV_long_covid_symptom_edge.csv": INCOV_PMIDS,
    "INSIGHT_OneFlorida_long_covid_symptom_edge.csv": INSIGHT_ONEFLORIDA_PMIDS,
    "RECOVER_long_covid_symptom_edge.csv": RECOVER_SYMPTOM_PMIDS,
}
PROTEIN_FILE = "RECOVER_PASC_protein_disease_edge.csv"
COMPOUND_FILE = "INCOV_long_covid_compound_edge.csv"

COL_DISEASE_ID = "disease_id"
COL_SYMPTOM_ID = "symptom_id"
COL_SYMPTOM_NAME = "symptom_name"
COL_PROTEIN_ID = "protein_id"
COL_PROTEIN_NAME = "protein_name"
COL_GENE_NAME = "gene_name"
COL_EDGE_TYPE = "edge_type"
COL_PMIDS = "pmid_list"
COL_PREPRINTS = "preprint_list"
COL_COMPOUND_ID = "compound_id"

# Columns that become node ids/names or edge publications; every OTHER column is kept verbatim as an edge attribute
SYMPTOM_NODE_COLS = {COL_DISEASE_ID, COL_SYMPTOM_ID, COL_SYMPTOM_NAME}
PROTEIN_NODE_COLS = {COL_DISEASE_ID, COL_PROTEIN_ID, COL_PROTEIN_NAME, COL_GENE_NAME, COL_PMIDS}
COMPOUND_NODE_COLS = {COL_DISEASE_ID, COL_COMPOUND_ID}

SOURCE_FILES_ATTR = "source_files"  # which of the tables back an edge
LIST_SPLIT = re.compile(r"[;,|\s]+")


class LongCovidHarmonizer(BaseHarmonizer):
    """Harmonizer for the ISB long COVID tables (Qi Wei, ISB; RECOVER-funded), originally imported into SPOKE.

    Five small CSVs, every row anchored on long COVID (DOID:0080848):
      * three symptom tables (INCOV, RECOVER, INSIGHT/OneFlorida) -> disease has_phenotype symptom (MeSH)
      * RECOVER protein table -> disease positively/negatively_correlated_with protein (UniProt), per edge_type
      * INCOV compound table -> compound applied_to_treat disease (ChEMBL)
    All edges are statistical_association / data_analysis_pipeline: dataset-specific cohort statistics.

    Every source column is retained: ids and names on nodes (gene symbol as a protein synonym), PMIDs as edge
    publications, and every other column verbatim as an edge attribute. A symptom reported by several tables is ONE
    edge carrying every table's columns (their names don't collide) plus the union of their PMIDs; `source_files`
    records which tables back it. Rows are kept exactly as given -- including MeSH terms that aren't really
    symptoms (e.g. Urinary Bladder, Sputum; typed NamedThing) and the compound table, whose percentages duplicate
    the first rows of the INCOV symptom table; both flagged to ISB.
    """

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
            raise ValueError(f"{self.source_name} requires input_file (the long-covid directory)")
        input_dir = Path(input_file)
        logging.info(f"Harmonizing {self.source_name}: {input_dir} -> {nodes_output}, {edges_output}")

        self._nodes: dict[str, dict[str, Any]] = {}  # raw curie -> {"category", "names", "synonyms"}
        self._edges: dict[tuple[str, str, str], dict[str, Any]] = {}  # (subject, predicate, object) -> accumulator

        for filename, pmids in SYMPTOM_FILES.items():
            for row in self._read_csv(input_dir / filename, required=(COL_DISEASE_ID, COL_SYMPTOM_ID)):
                disease = self._add_node(row[COL_DISEASE_ID], DISEASE_CATEGORY)
                symptom_id = row[COL_SYMPTOM_ID]
                category = ROOT_CATEGORY if symptom_id in NON_PHENOTYPE_SYMPTOM_IDS else SYMPTOM_CATEGORY
                symptom = self._add_node(f"{SYMPTOM_VOCAB}:{symptom_id}", category, row.get(COL_SYMPTOM_NAME))
                self._add_edge(disease, HAS_PHENOTYPE, symptom, filename, row, SYMPTOM_NODE_COLS, pmids)

        n_bad_edge_type = 0
        for row in self._read_csv(input_dir / PROTEIN_FILE, required=(COL_DISEASE_ID, COL_PROTEIN_ID)):
            predicate = EDGE_TYPE_TO_PREDICATE.get(row.get(COL_EDGE_TYPE, ""))
            if predicate is None:
                n_bad_edge_type += 1
                logging.warning(f"{self.source_name}: unknown edge_type {row.get(COL_EDGE_TYPE)!r}; skipping row")
                continue
            disease = self._add_node(row[COL_DISEASE_ID], DISEASE_CATEGORY)
            protein = self._add_node(
                f"{PROTEIN_VOCAB}:{row[COL_PROTEIN_ID]}",
                PROTEIN_CATEGORY,
                row.get(COL_PROTEIN_NAME),
                synonym=row.get(COL_GENE_NAME),
            )
            pmids = [f"PMID:{pmid}" for pmid in self._split_list(row.get(COL_PMIDS))]
            self._add_edge(disease, predicate, protein, PROTEIN_FILE, row, PROTEIN_NODE_COLS, pmids)

        for row in self._read_csv(input_dir / COMPOUND_FILE, required=(COL_DISEASE_ID, COL_COMPOUND_ID)):
            disease = self._add_node(row[COL_DISEASE_ID], DISEASE_CATEGORY)
            compound = self._add_node(f"{COMPOUND_VOCAB}:{row[COL_COMPOUND_ID]}", COMPOUND_CATEGORY)
            self._add_edge(compound, APPLIED_TO_TREAT, disease, COMPOUND_FILE, row, COMPOUND_NODE_COLS, INCOV_PMIDS)

        nodes = [self._build_node(curie, info) for curie, info in self._nodes.items()]
        edges = [self._build_edge(key, acc) for key, acc in self._edges.items()]
        save_to_jsonl(nodes, nodes_output, mode="w")
        save_to_jsonl(edges, edges_output, mode="w")
        logging.info(
            f"{self.source_name} harmonization complete: {len(nodes)} nodes, {len(edges)} edges "
            f"({n_bad_edge_type} protein rows skipped for unknown edge_type)"
        )

    # ------------------------------------------------------------------ readers

    @staticmethod
    def _read_csv(path: Path, required: tuple[str, ...]):
        """Rows with every `required` column filled, values stripped. Skips the trailing blank rows the tables
        end with."""
        with open(path, encoding="utf-8-sig", newline="") as f:
            for row in csv.DictReader(f):
                row = {key: (value or "").strip() for key, value in row.items() if key}
                if all(row.get(col) for col in required):
                    yield row

    # ------------------------------------------------------------------ accumulators

    def _add_node(self, curie: str, category: str, name: str | None = None, synonym: str | None = None) -> str:
        """Names are ordered: the first becomes the node's name, the rest (and any synonym) its synonyms."""
        info = self._nodes.setdefault(curie, {"category": category, "names": [], "synonyms": []})
        if name and name not in info["names"]:
            info["names"].append(name)
        if synonym and synonym not in info["synonyms"]:
            info["synonyms"].append(synonym)
        return curie

    def _add_edge(
        self,
        subject: str,
        predicate: str,
        obj: str,
        filename: str,
        row: dict,
        node_cols: set[str],
        publications: list[str],
    ):
        acc = self._edges.setdefault(
            (subject, predicate, obj), {"attributes": {SOURCE_FILES_ATTR: []}, "publications": []}
        )
        attributes = acc["attributes"]
        if filename not in attributes[SOURCE_FILES_ATTR]:
            attributes[SOURCE_FILES_ATTR].append(filename)
        for col, raw in row.items():
            if col in node_cols or raw == "":
                continue
            value = self._parse_value(col, raw)
            existing = attributes.get(col)
            if existing is None:
                attributes[col] = value
            elif existing != value:
                # The same edge reported twice in one table (e.g. a protein in two papers): keep both values
                merged = existing if isinstance(existing, list) else [existing]
                for item in value if isinstance(value, list) else [value]:
                    if item not in merged:
                        merged.append(item)
                attributes[col] = merged
        for pmid in publications:
            if pmid not in acc["publications"]:
                acc["publications"].append(pmid)

    # ------------------------------------------------------------------ builders

    def _build_node(self, curie: str, info: dict) -> dict:
        names = info["names"]
        synonyms = list(dict.fromkeys(names[1:] + info["synonyms"]))
        return self.create_node(
            curie=curie,
            categories=[info["category"]],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=names[0] if names else None,
            synonyms=synonyms or None,
        )

    def _build_edge(self, key: tuple[str, str, str], acc: dict) -> dict:
        subject, predicate, obj = key
        return self.create_edge(
            subject_id=subject,
            object_id=obj,
            predicate=predicate,
            primary_ks=self.source_infores,
            knowledge_level=STATISTICAL_ASSOCIATION,
            agent_type=DATA_ANALYSIS_PIPELINE,
            publications=acc["publications"] or None,
            attributes=acc["attributes"],
        )

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _split_list(value: str | None) -> list[str]:
        return [item for item in LIST_SPLIT.split(value or "") if item]

    @classmethod
    def _parse_value(cls, col: str, raw: str):
        if col == COL_PREPRINTS:
            return cls._split_list(raw)
        try:
            return float(raw)
        except ValueError:
            return raw
