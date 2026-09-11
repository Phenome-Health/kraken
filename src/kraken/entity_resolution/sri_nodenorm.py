"""SRI Node Normalizer client for names + categories on bare CURIEs (plan §4).

Bare ids arriving in an equivalency list have no name/category from their
source; without a category nearly every leaf node is a wildcard and the branch
guardrail is inert. This client supplies both, using the RENCI API directly.

Key decisions (correcting the plan doc, per Amy):

* ``conflate=true`` (gene/protein conflation) and ``drug_chemical_conflate=true``
  (drug/chemical conflation) — both wanted. NOTE: the cache key is the CURIE only,
  so changing either flag makes prior cache entries stale — clear the cache and
  re-resolve after a flag change.
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
import time
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
    "MONDO": "biolink:Disease",
    # HP terms are usually phenotypes but sometimes diseases, so back off to the shared parent.
    "HP": "biolink:DiseaseOrPhenotypicFeature",
    "UBERON": "biolink:AnatomicalEntity",
    "CL": "biolink:Cell",
    # NOTE: DRUGBANK is intentionally NOT here — it spans small molecules AND biologics/
    # protein drugs, so no single category is safe (and this last-resort backup is rarely hit).
    # NOTE: REACT is intentionally NOT here. Reactome R-HSA ids are ambiguous
    # (pathways, reactions, AND physical entities like modified-protein states), so a
    # blanket REACT->Pathway guess mis-types many nodes (e.g. Reactome TP53 PTM-forms
    # that KG2 lists as gene equivalents). Same reason ENSEMBL / bare GO / KEGG are omitted.
}

# Prefixes the normalizer cannot resolve (structural strings, not identifiers):
# skip them entirely rather than wasting an API round trip.
NON_QUERYABLE_PREFIXES: frozenset[str] = frozenset({"SMILES", "INCHI"})

# Prefix -> taxon for single-species nomenclature authorities (the prefix itself
# DEFINES the species). Used ONLY as a backup when the normalizer returns no
# taxon (e.g. an id the normalizer doesn't recognize) — never overriding a real
# NN answer. Kept to prefixes that are unambiguously one species; multi-species
# prefixes (NCBIGene, UniProtKB, ENSEMBL, Xenbase) are deliberately omitted.
TAXON_BY_PREFIX: dict[str, str] = {
    "HGNC": "NCBITaxon:9606",       # human
    "MGI": "NCBITaxon:10090",       # mouse
    "RGD": "NCBITaxon:10116",       # rat
    "ZFIN": "NCBITaxon:7955",       # zebrafish
    "FB": "NCBITaxon:7227",         # fruit fly (FlyBase)
    "FlyBase": "NCBITaxon:7227",
    "WB": "NCBITaxon:6239",         # C. elegans (WormBase)
    "WormBase": "NCBITaxon:6239",
    "SGD": "NCBITaxon:559292",      # S. cerevisiae S288C
    "PomBase": "NCBITaxon:4896",    # S. pombe
    "dictyBase": "NCBITaxon:44689",  # D. discoideum
    "TAIR": "NCBITaxon:3702",       # A. thaliana
}


@dataclass(frozen=True)
class NormInfo:
    """Per-CURIE facts from the normalizer. ``taxa`` and ``categories`` may be
    empty; ``canonical`` is the clique's canonical id (used to group cliques for
    equivalence evidence). A member shares a clique with every other member that
    has the same ``canonical``."""

    label: str | None
    categories: tuple[str, ...]
    taxa: tuple[str, ...] = ()
    canonical: str | None = None


def infer_category(curie: str) -> str | None:
    return PREFIX_CATEGORY.get(curie.split(":", 1)[0])


def infer_taxon(curie: str) -> str | None:
    """Backup taxon from a single-species nomenclature prefix (or None)."""
    return TAXON_BY_PREFIX.get(curie.split(":", 1)[0])


class NodeNormClient:
    def __init__(
        self,
        cache_path: str | Path,
        *,
        base_url: str = DEFAULT_BASE_URL,
        conflate: bool = True,
        drug_chemical_conflate: bool = True,
        individual_types: bool = True,
        batch_size: int = 1000,
        timeout: float = 60.0,
        max_retries: int = 4,
        retry_backoff: float = 2.0,
        progress_every: int = 200,
        session: requests.Session | None = None,
    ):
        self.base_url = base_url
        self.conflate = conflate
        self.drug_chemical_conflate = drug_chemical_conflate
        self.individual_types = individual_types
        self.batch_size = batch_size
        self.timeout = timeout
        self.max_retries = max_retries  # retry a failed batch this many times (exponential backoff)
        self.retry_backoff = retry_backoff  # base seconds; sleep = retry_backoff * 2**attempt
        self.progress_every = progress_every  # log + commit cache every N batches
        self._session = session or requests.Session()
        self._db = sqlite3.connect(str(cache_path))
        self._db.execute(
            "CREATE TABLE IF NOT EXISTS norm_cache ("
            "curie TEXT PRIMARY KEY, canonical TEXT, label TEXT, categories TEXT, "
            "taxa TEXT, resolved INTEGER)"
        )
        # Migrate an older cache (label/categories/resolved only) by adding the new
        # columns; rows written before this change lack canonical/taxa, so a rerun
        # that needs cliques or taxa must be run against a cleared cache.
        existing = {r[1] for r in self._db.execute("PRAGMA table_info(norm_cache)")}
        for column in ("canonical", "taxa"):
            if column not in existing:
                self._db.execute(f"ALTER TABLE norm_cache ADD COLUMN {column} TEXT")
        self._db.execute("CREATE INDEX IF NOT EXISTS idx_norm_canonical ON norm_cache(canonical)")
        self._db.commit()

    # ---- cache ----

    def _cache_get(self, curie: str) -> NormInfo | None:
        row = self._db.execute(
            "SELECT label, categories, taxa, canonical FROM norm_cache WHERE curie = ?", (curie,)
        ).fetchone()
        if row is None:
            return None
        label, categories_json, taxa_json, canonical = row
        categories = tuple(json.loads(categories_json)) if categories_json else ()
        taxa = tuple(json.loads(taxa_json)) if taxa_json else ()
        return NormInfo(label=label, categories=categories, taxa=taxa, canonical=canonical)

    def _cache_put(self, curie: str, info: NormInfo, resolved: bool) -> None:
        self._db.execute(
            "INSERT OR REPLACE INTO norm_cache (curie, canonical, label, categories, taxa, resolved) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            (
                curie,
                info.canonical,
                info.label,
                json.dumps(list(info.categories)),
                json.dumps(list(info.taxa)),
                int(resolved),
            ),
        )

    # ---- API ----

    def _params(self) -> dict:
        return {
            "conflate": str(self.conflate).lower(),
            "drug_chemical_conflate": str(self.drug_chemical_conflate).lower(),
            "individual_types": str(self.individual_types).lower(),
        }

    def _fetch_batch(self, curies: list[str]) -> dict[str, NormInfo] | None:
        """Fetch one batch, retrying transient failures with exponential backoff.
        Returns the parsed dict on success (curies absent from it are legitimately
        unresolved), or ``None`` if it failed after all retries (caller must NOT
        cache those — they should be retried on the next run, not poisoned empty)."""
        payload = {"curies": curies, **{k: v == "true" for k, v in self._params().items()}}
        for attempt in range(self.max_retries + 1):
            try:
                resp = self._session.post(self.base_url, json=payload, timeout=self.timeout)
                resp.raise_for_status()
                return self._parse_response(resp.json())
            except (requests.RequestException, ValueError) as exc:
                if attempt < self.max_retries:
                    sleep_s = self.retry_backoff * (2**attempt)
                    logging.warning(
                        "Node Normalizer batch failed (attempt %d/%d), retrying in %.0fs: %s",
                        attempt + 1,
                        self.max_retries + 1,
                        sleep_s,
                        exc,
                    )
                    time.sleep(sleep_s)
                else:
                    logging.error(
                        "Node Normalizer batch of %d curies failed after %d attempts; leaving uncached "
                        "(will retry next run): %s",
                        len(curies),
                        self.max_retries + 1,
                        exc,
                    )
        return None

    @staticmethod
    def _parse_response(data: dict) -> dict[str, NormInfo]:
        """HARVEST every clique member from the response, not just the queried ids.

        One query returns the full clique (every ``equivalent_identifiers`` row
        with its own label/type/taxa), so we record all of them — the caller
        caches them, which both de-duplicates future queries (clique-mates are
        already resolved) and provides the clique membership (via ``canonical``)
        used to emit equivalence evidence. Per-id ``type``/``taxa`` fall back to
        the clique-level values when a row omits them."""
        harvested: dict[str, NormInfo] = {}
        for entry in data.values():
            if not entry:
                continue
            canonical = (entry.get("id") or {}).get("identifier")
            clique_types = entry.get("type") or []
            clique_taxa = entry.get("taxa") or []
            for eq in entry.get("equivalent_identifiers", []):
                member = eq.get("identifier")
                if not member:
                    continue
                types = eq.get("type") or clique_types
                categories = tuple(types) if isinstance(types, list) else (types,)
                taxa_raw = eq.get("taxa") or clique_taxa
                taxa = tuple(taxa_raw) if isinstance(taxa_raw, list) else (taxa_raw,)
                harvested[member] = NormInfo(
                    label=eq.get("label"),
                    categories=categories,
                    taxa=taxa,
                    canonical=canonical,
                )
        return harvested

    # ---- public ----

    def resolve(self, curies: Iterable[str]) -> dict[str, NormInfo]:
        """Return ``{curie: NormInfo}`` holding ONLY the normalizer's own facts (categories,
        taxa, label, clique). Prefix backups are NOT applied here: a category/taxon prefix
        guess is a last resort that must never pre-empt a source-derived value, so the
        CALLER (entity_resolution/build.py) applies ``infer_category`` / ``infer_taxon``
        only AFTER source values. Non-queryable prefixes (SMILES etc.) skip the API and
        resolve to empty facts.

        HARVEST + DEDUP: one query returns a CURIE's whole clique, so every
        clique-mate's facts are cached too. A queued id already filled in by an
        earlier batch's harvest is skipped, so each clique is fetched once, not
        once per member. The cached ``canonical`` also lets ``iter_cliques`` emit
        equivalence evidence from the normalizer's cliques.
        """
        wanted = list(dict.fromkeys(curies))  # dedup input, keep order
        result: dict[str, NormInfo] = {}
        queue: list[str] = []
        cached_hits = non_queryable = 0

        for curie in wanted:
            cached = self._cache_get(curie)
            if cached is not None:
                result[curie] = cached
                cached_hits += 1
            elif curie.split(":", 1)[0] in NON_QUERYABLE_PREFIXES:
                # SMILES/INCHI are raw structure strings the normalizer cannot resolve
                # (verified) -- skip the API; they resolve to empty facts (any prefix backup
                # is applied by the caller, after source values).
                info = NormInfo(None, ())
                result[curie] = info
                self._cache_put(curie, info, resolved=False)
                non_queryable += 1
            else:
                queue.append(curie)

        logging.info(
            "Node Normalizer: up to %d curies to fetch (%d already cached, %d non-queryable "
            "[SMILES/INCHI, skipped]), harvesting cliques",
            len(queue),
            cached_hits,
            non_queryable,
        )
        pending: list[str] = []
        done = fetched_batches = failed_batches = 0

        def flush_batch() -> None:
            nonlocal done, fetched_batches, failed_batches
            if not pending:
                return
            fetched_batches += 1
            harvested = self._fetch_batch(pending)
            if harvested is None:
                failed_batches += 1
                for c in pending:  # transient failure this run: empty facts, do NOT cache
                    result[c] = NormInfo(None, ())
            else:
                # cache EVERY harvested member (clique-mates included) so future
                # queued ids in the same clique become cache hits (the dedup). We store
                # NN's own facts only -- prefix backups are the caller's last resort.
                for member, info in harvested.items():
                    self._cache_put(member, info, resolved=True)
                for c in pending:  # a queried id NN didn't return -> unrecognized
                    result[c] = self._cache_get(c) or self._record_unrecognized(c)
            done += len(pending)
            if fetched_batches % self.progress_every == 0:
                self._db.commit()
                logging.info(
                    "Node Normalizer: fetched %d curies (%d batches%s)",
                    done,
                    fetched_batches,
                    f", {failed_batches} failed" if failed_batches else "",
                )
            pending.clear()

        for curie in queue:
            cached = self._cache_get(curie)  # may have been harvested by an earlier batch
            if cached is not None:
                result[curie] = cached
                continue
            pending.append(curie)
            if len(pending) >= self.batch_size:
                flush_batch()
        flush_batch()

        self._db.commit()
        return result

    def get(self, curie: str) -> NormInfo | None:
        """Cached facts for one CURIE (or None if never resolved). Public accessor
        so materialization can look up an isolated id's category/label/taxa."""
        return self._cache_get(curie)

    def iter_labels(self):
        """Yield ``(curie, label)`` for every cached id that has a normalizer label.
        This is the PER-ID name source for name-similarity: every identifier is named
        individually (by its own NN label), not by some source node's primary name."""
        yield from self._db.execute("SELECT curie, label FROM norm_cache WHERE label IS NOT NULL AND label != ''")

    def iter_cliques(self):
        """Yield ``(canonical, [member curies])`` for every normalizer clique in the
        cache (resolved rows only), streamed in canonical order so memory is bounded
        by one clique. This is the equivalence signal for the match graph."""
        cursor = self._db.execute(
            "SELECT canonical, curie FROM norm_cache "
            "WHERE canonical IS NOT NULL AND resolved = 1 ORDER BY canonical"
        )
        current: str | None = None
        members: list[str] = []
        for canonical, curie in cursor:
            if canonical != current and members:
                yield current, members
                members = []
            current = canonical
            members.append(curie)
        if members:
            yield current, members

    def _record_unrecognized(self, curie: str) -> NormInfo:
        """Cache + return empty facts for an id the normalizer couldn't resolve.

        NO prefix backup is applied here. Category/taxon prefix guesses are a LAST resort
        that must never pre-empt a source-derived value, so they are applied by the CALLER
        (entity_resolution/build.py) only AFTER source values — via ``infer_category`` /
        ``infer_taxon``. This keeps the normalizer client a pure record of what NN knows."""
        info = NormInfo(None, ())
        self._cache_put(curie, info, resolved=False)
        return info

    def close(self) -> None:
        self._db.commit()
        self._db.close()

    def __enter__(self) -> NodeNormClient:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()
