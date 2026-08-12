# pgs_catalog.py
import logging
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import openpyxl

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import COMPUTATIONAL_MODEL, PGS_CATALOG_SOURCE_ID, STATISTICAL_ASSOCIATION
from kraken.utils.kg_io import save_to_jsonl

# --- v1 selection knobs ---
# Rather than a high global evaluation threshold (which silently drops important-but-less-studied traits like
# asthma/Alzheimer), we keep the best-validated score PER TRAIT. This keeps every trait's flagship PGS while
# cutting redundancy (e.g. ~22 CAD scores -> the single most-evaluated one).
MIN_EVALUATION_SAMPLE_SETS = 5  # a PGS must be validated in >= this many distinct evaluation sample sets
TOP_N_PER_TRAIT = 1  # keep only the N best-validated (most-evaluated) scores per mapped trait

# --- PLACEHOLDER Biolink types (pending a types review; each is a single swappable constant) ---
# There is no clean Biolink class for a polygenic score; InformationContentEntity is a stand-in.
PGS_NODE_CATEGORY = "biolink:InformationContentEntity"  # TODO(types): revisit (ClinicalFinding? a custom PGS class?)
TRAIT_STUB_CATEGORY = "biolink:NamedThing"  # minimal stub; real category arrives when ontology sources merge in
PGS_TRAIT_PREDICATE = "biolink:related_to"  # TODO(types): revisit (e.g. assesses_risk_of / predicts)
# NOTE: PGS -> gene edges (positional annotation of the scoring files' variants) are a planned later increment;
# see the TODO in harmonize(). The scoring-file FTP link is retained on each PGS node's attributes for that step.

PGS_CURIE_PREFIX = "PGS"  # PGS Catalog accession is e.g. "PGS000013"; emitted as "PGS:000013"


class PGSCatalogHarmonizer(BaseHarmonizer):
    """Harmonizer for the PGS Catalog metadata bundle (``pgs_all_metadata.xlsx``).

    v1 scope: from the ~7k scores, keep the single best-validated PGS per mapped trait (ranked by number of
    distinct evaluation sample sets, gated at >= ``MIN_EVALUATION_SAMPLE_SETS``). Each selected score becomes a
    PGS node linked to its mapped EFO/MONDO trait(s); a minimal stub node is minted for each trait endpoint so
    nothing orphans. PGS -> gene edges (via positional annotation of the per-PGS scoring files' variants) are a
    planned later increment and are NOT built here.

    All node/edge Biolink types are placeholders -- see the module-level ``*_CATEGORY`` / ``*_PREDICATE``
    constants (pending a types review); each is a single swappable constant.
    """

    source_infores = PGS_CATALOG_SOURCE_ID

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

        wb = openpyxl.load_workbook(input_file, read_only=True)
        scores = self._load_scores(wb)  # pgs_id -> metadata dict
        evaluations = self._load_evaluations(wb)  # pgs_id -> set of distinct evaluation sample sets
        wb.close()

        selected = self._select_top_per_trait(scores, evaluations)

        pgs_nodes: list[dict] = []
        edges: list[dict] = []
        trait_stubs: dict[str, dict] = {}
        for pgs_id in selected:
            meta = scores[pgs_id]
            pgs_nodes.append(self._make_pgs_node(pgs_id, meta))
            for trait_curie, trait_label in meta["trait_curies"]:
                trait_stubs.setdefault(trait_curie, self._make_trait_stub(trait_curie, trait_label))
                edges.append(self._make_pgs_trait_edge(pgs_id, trait_curie, len(evaluations.get(pgs_id, ()))))

        all_nodes = pgs_nodes + list(trait_stubs.values())
        save_to_jsonl(all_nodes, nodes_output, mode="w")
        save_to_jsonl(edges, edges_output, mode="w")
        logging.info(
            f"{self.source_name} harmonization complete: {len(pgs_nodes)} PGS nodes + {len(trait_stubs)} trait "
            f"stubs, {len(edges)} PGS->trait edges (selected {len(selected)} of {len(scores)} scores)"
        )

    # ------------------------------------------------------------------ loading

    def _load_scores(self, wb) -> dict[str, dict[str, Any]]:
        scores: dict[str, dict] = {}
        for r in wb["Scores"].iter_rows(min_row=2, values_only=True):
            pgs_id = r[0]
            if not pgs_id:
                continue
            scores[pgs_id] = {
                "name": r[1],
                "reported_trait": r[2],
                "mapped_trait_label": r[3],
                "mapped_trait_efo": r[4],
                "trait_curies": self._parse_trait_curies(r[3], r[4]),
                "development_method": r[5],
                "genome_build": r[7],
                "num_variants": self._to_int(r[8]),
                "weight_type": r[10],
                "pgp_id": r[11],
                "pmid": r[12],
                "ftp_link": r[18],  # scoring file -- retained for the later gene-annotation increment
                "release_date": str(r[19]) if r[19] else None,  # str(): cell may be a datetime (not JSON-safe)
            }
        return scores

    def _load_evaluations(self, wb) -> dict[str, set]:
        """pgs_id -> set of distinct evaluation sample sets (PSS). Distinct sample sets, not raw metric rows,
        is the honest 'independently validated in N cohorts' measure (one cohort often reports several metrics)."""
        sample_sets: dict[str, set] = defaultdict(set)
        for r in wb["Performance Metrics"].iter_rows(min_row=2, values_only=True):
            pgs_id, sample_set = r[1], r[2]
            if pgs_id and sample_set:
                sample_sets[pgs_id].add(sample_set)
        return sample_sets

    # ------------------------------------------------------------------ selection

    def _select_top_per_trait(self, scores: dict, evaluations: dict[str, set]) -> set[str]:
        eval_count = lambda pgs_id: len(evaluations.get(pgs_id, ()))
        # Group eligible scores by each mapped trait, then keep the top-N most-evaluated per trait.
        by_trait: dict[str, list[str]] = defaultdict(list)
        for pgs_id, meta in scores.items():
            if eval_count(pgs_id) >= MIN_EVALUATION_SAMPLE_SETS:
                for trait_curie, _ in meta["trait_curies"]:
                    by_trait[trait_curie].append(pgs_id)
        selected: set[str] = set()
        for pgs_ids in by_trait.values():
            selected.update(sorted(pgs_ids, key=eval_count, reverse=True)[:TOP_N_PER_TRAIT])
        return selected

    # ------------------------------------------------------------------ node/edge builders

    def _make_pgs_node(self, pgs_id: str, meta: dict) -> dict:
        curie = self._pgs_curie(pgs_id)
        publications = [f"PMID:{meta['pmid']}"] if meta.get("pmid") else None
        attributes = {
            "pgs_catalog_id": pgs_id,  # native accession, e.g. "PGS000013"
            "reported_trait": meta["reported_trait"],
            "mapped_trait_label": meta["mapped_trait_label"],
            "mapped_trait_efo": meta["mapped_trait_efo"],
            "pgs_development_method": meta["development_method"],
            "genome_build": meta["genome_build"],
            "number_of_variants": meta["num_variants"],
            "type_of_variant_weight": meta["weight_type"],
            "pgs_publication_id": meta["pgp_id"],
            "scoring_file_ftp": meta["ftp_link"],
            "release_date": meta["release_date"],
        }
        return self.create_node(
            curie=curie,
            categories=[PGS_NODE_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=[curie],
            name=meta["name"],
            description=f"Polygenic score for {meta['reported_trait']}" if meta.get("reported_trait") else None,
            publications=publications,
            attributes={k: v for k, v in attributes.items() if v not in (None, "")},
        )

    def _make_trait_stub(self, trait_curie: str, trait_label: str | None) -> dict:
        return self.create_node(
            curie=trait_curie,
            categories=[TRAIT_STUB_CATEGORY],
            provided_by=self.source_infores,
            equivalent_ids=[trait_curie],
            name=trait_label,
        )

    def _make_pgs_trait_edge(self, pgs_id: str, trait_curie: str, num_sample_sets: int) -> dict:
        return self.create_edge(
            subject_id=self._pgs_curie(pgs_id),
            object_id=trait_curie,
            predicate=PGS_TRAIT_PREDICATE,
            primary_ks=self.source_infores,
            knowledge_level=STATISTICAL_ASSOCIATION,
            agent_type=COMPUTATIONAL_MODEL,
            attributes={"num_evaluation_sample_sets": num_sample_sets} if num_sample_sets else None,
        )

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _pgs_curie(pgs_id: str) -> str:
        return f"{PGS_CURIE_PREFIX}:{pgs_id.removeprefix('PGS')}"  # PGS000013 -> PGS:000013

    @staticmethod
    def _parse_trait_curies(label: str | None, efo: str | None) -> list[tuple[str, str | None]]:
        """Map the mapped-trait EFO id(s) to (curie, label) pairs, e.g. 'MONDO_0005010' -> ('MONDO:0005010', ...).
        Both columns may hold multiple '|'/','-separated values; labels are zipped positionally when counts line up."""
        if not efo:
            return []
        efo_ids = [e.strip() for e in re.split(r"[|,]", str(efo)) if e.strip()]
        labels = [l.strip() for l in re.split(r"[|,]", str(label))] if label else []
        pairs = []
        for i, raw in enumerate(efo_ids):
            curie = raw.replace("_", ":", 1)  # MONDO_0005010 -> MONDO:0005010
            pairs.append((curie, labels[i] if i < len(labels) else (label or None)))
        return pairs

    @staticmethod
    def _to_int(value: Any) -> int | None:
        try:
            return int(str(value).replace(",", "")) if value not in (None, "") else None
        except (ValueError, TypeError):
            return None
