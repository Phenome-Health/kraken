"""SRI Node Normalizer client for names + categories on bare CURIEs (plan §4).

Bare ids arriving in an equivalency list have no name/category from their
source; without a category nearly every leaf node is a wildcard and the branch
guardrail is inert. This client supplies both, using the RENCI API directly.

Key decisions (correcting the plan doc, per Amy):

* ``conflate=true`` — we WANT gene/protein conflation — and
  ``drug_chemical_conflate=false``.
* ``individual_types=true`` so each equivalent identifier carries its own
  category (not just the clique's).
* Per-CURIE label comes from that CURIE's entry in ``equivalent_identifiers``,
  **not** the top-level ``id.label`` (which is the clique's preferred label and
  would smuggle Babel's clustering back in via naming).
* **The normalizer is the source of truth for categories.** Prefix inference
  (HGNC->Gene, ...) is only a **backup** for ids the normalizer doesn't recognize
  — a prefix->category heuristic must never override a real answer.
* **Persistent on-disk cache** keyed by CURIE (sqlite), including negatives, so
  the (higher, now-query-everything) request volume is paid once.
* Batches of 1000.
"""

from __future__ import annotations

import json
import logging
import sqlite3
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import requests

DEFAULT_BASE_URL = "https://nodenormalization-sri.renci.org/get_normalized_nodes"

# Unambiguous prefix -> Biolink category. Conservative on purpose (ENSEMBL and
# bare GO/KEGG are ambiguous, so omitted). Extend as safe. Inferable prefixes
# skip the API entirely.
PREFIX_CATEGORY: dict[str, str] = {
    "HGNC": "biolink:Gene",
    "NCBIGene": "biolink:Gene",
    "UniProtKB": "biolink:Protein",
    "PR": "biolink:Protein",
    "NCBITaxon": "biolink:OrganismTaxon",
    "CHEBI": "biolink:ChemicalEntity",
    "PUBCHEM.COMPOUND": "biolink:SmallMolecule",
    "KEGG.COMPOUND": "biolink:SmallMolecule",
    "HMDB": "biolink:SmallMolecule",
    "INCHIKEY": "biolink:ChemicalEntity",
    "DRUGBANK": "biolink:Drug",
    "MONDO": "biolink:Disease",
    "HP": "biolink:PhenotypicFeature",
    "UBERON": "biolink:AnatomicalEntity",
    "CL": "biolink:Cell",
    "REACT": "biolink:Pathway",
}

# Prefixes the normalizer cannot resolve (structural strings, not identifiers):
# skip them entirely rather than wasting an API round trip.
NON_QUERYABLE_PREFIXES: frozenset[str] = frozenset({"SMILES", "INCHI"})


@dataclass(frozen=True)
class NormInfo:
    """Label + categories for one CURIE (categories may be empty)."""

    label: str | None
    categories: tuple[str, ...]


def infer_category(curie: str) -> str | None:
    return PREFIX_CATEGORY.get(curie.split(":", 1)[0])


class NodeNormClient:
    def __init__(
        self,
        cache_path: str | Path,
        *,
        base_url: str = DEFAULT_BASE_URL,
        conflate: bool = True,
        drug_chemical_conflate: bool = False,
        individual_types: bool = True,
        batch_size: int = 1000,
        timeout: float = 60.0,
        session: requests.Session | None = None,
    ):
        self.base_url = base_url
        self.conflate = conflate
        self.drug_chemical_conflate = drug_chemical_conflate
        self.individual_types = individual_types
        self.batch_size = batch_size
        self.timeout = timeout
        self._session = session or requests.Session()
        self._db = sqlite3.connect(str(cache_path))
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS norm_cache ("
            "curie TEXT PRIMARY KEY, label TEXT, categories TEXT, resolved INTEGER)"
        )
        self._db.commit()

    # ---- cache ----

    def _cache_get(self, curie: str) -> NormInfo | None:
        row = self._db.execute(
            "SELECT label, categories, resolved FROM norm_cache WHERE curie = ?", (curie,)
        ).fetchone()
        if row is None:
            return None
        label, categories_json, _resolved = row
        categories = tuple(json.loads(categories_json)) if categories_json else ()
        return NormInfo(label=label, categories=categories)

    def _cache_put(self, curie: str, info: NormInfo, resolved: bool) -> None:
        self._db.execute(
            "INSERT OR REPLACE INTO norm_cache (curie, label, categories, resolved) VALUES (?, ?, ?, ?)",
            (curie, info.label, json.dumps(list(info.categories)), int(resolved)),
        )

    # ---- API ----

    def _params(self) -> dict:
        return {
            "conflate": str(self.conflate).lower(),
            "drug_chemical_conflate": str(self.drug_chemical_conflate).lower(),
            "individual_types": str(self.individual_types).lower(),
        }

    def _fetch_batch(self, curies: list[str]) -> dict[str, NormInfo]:
        payload = {"curies": curies, **{k: v == "true" for k, v in self._params().items()}}
        try:
            resp = self._session.post(self.base_url, json=payload, timeout=self.timeout)
            resp.raise_for_status()
            data = resp.json()
        except (requests.RequestException, ValueError) as exc:
            logging.warning("Node Normalizer request failed for %d curies: %s", len(curies), exc)
            return {}
        return self._parse_response(curies, data)

    @staticmethod
    def _parse_response(curies: Iterable[str], data: dict) -> dict[str, NormInfo]:
        out: dict[str, NormInfo] = {}
        for curie in curies:
            entry = data.get(curie)
            if not entry:
                continue
            # per-CURIE label + type from the matching equivalent_identifiers row
            label: str | None = None
            categories: tuple[str, ...] = ()
            for eq in entry.get("equivalent_identifiers", []):
                if eq.get("identifier") == curie:
                    label = eq.get("label")
                    types = eq.get("type") or []
                    categories = tuple(types) if isinstance(types, list) else (types,)
                    break
            if not categories:
                # fall back to the clique-level type only for the category
                top_types = entry.get("type") or []
                categories = tuple(top_types) if isinstance(top_types, list) else (top_types,)
            out[curie] = NormInfo(label=label, categories=categories)
        return out

    # ---- public ----

    def resolve(self, curies: Iterable[str], *, use_inference: bool = True) -> dict[str, NormInfo]:
        """Return ``{curie: NormInfo}``. The normalizer is the source of truth:
        every id is queried (batched) and its result cached, including negatives.
        Prefix inference is used ONLY as a **backup category** when the normalizer
        doesn't recognize an id (or returns no type) — trusting a prefix→category
        heuristic over the normalizer is exactly the kind of shortcut that goes
        wrong, so it never overrides a real answer. Non-queryable prefixes (raw
        structure strings like SMILES) skip the API and fall straight to inference.
        """
        wanted = list(dict.fromkeys(curies))  # dedup, keep order
        result: dict[str, NormInfo] = {}
        to_fetch: list[str] = []

        for curie in wanted:
            cached = self._cache_get(curie)
            if cached is not None:
                result[curie] = cached
                continue
            if curie.split(":", 1)[0] in NON_QUERYABLE_PREFIXES:
                info = self._inference_backup(curie) if use_inference else NormInfo(label=None, categories=())
                result[curie] = info
                self._cache_put(curie, info, resolved=False)
                continue
            to_fetch.append(curie)

        for start in range(0, len(to_fetch), self.batch_size):
            batch = to_fetch[start : start + self.batch_size]
            fetched = self._fetch_batch(batch)
            for curie in batch:
                info = fetched.get(curie)
                recognized = info is not None and bool(info.categories)
                if not recognized and use_inference:
                    # normalizer didn't type it -> prefix inference as a backup,
                    # keeping any label the normalizer did return.
                    backup = self._inference_backup(curie)
                    label = info.label if info is not None else None
                    info = NormInfo(label=label, categories=backup.categories)
                elif info is None:
                    info = NormInfo(label=None, categories=())
                result[curie] = info
                self._cache_put(curie, info, resolved=recognized)

        self._db.commit()
        return result

    @staticmethod
    def _inference_backup(curie: str) -> NormInfo:
        inferred = infer_category(curie)
        return NormInfo(label=None, categories=(inferred,) if inferred else ())

    def close(self) -> None:
        self._db.commit()
        self._db.close()

    def __enter__(self) -> NodeNormClient:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
