"""Materialize a cluster into one canonical node (plan §1, §7 representative).

Order-independence is required (the whole point of the overhaul), so every
reconciliation here is a function of the *set* of members, not their order:

* union fields (equivalent_ids, synonyms, provided_by, urls, publications) are
  set-unions;
* the representative id / name / description are chosen by a **total order**
  (ranked prefix for the cluster's family, then CURIE), which is
  order-independent; a few chemical names are corrected by hand (NAME_CORRECTIONS);
* categories are the leaf-filtered union, plus Drug for a chemical holding an ATC drug code;
* taxon is the single guardrail-permitted taxon (or None on conflict).

The representative id is cosmetic (queries resolve through the membership map),
so it only labels the cluster; name and description may come from different
members (a lower-ranked id can supply the description the top id lacks).
"""

from __future__ import annotations

import logging
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import yaml

from kraken.entity_resolution.families import ALL_FAMILIES, CHEMICAL_FAMILY, BranchFamilies
from kraken.utils.constants import (
    NODE_ATTRIBUTES,
    NODE_CATEGORIES,
    NODE_CHEMICAL_FORMULA,
    NODE_DESCRIPTION,
    NODE_EQUIVALENT_IDS,
    NODE_EXACT_MASS,
    NODE_ID,
    NODE_NAME,
    NODE_PROVIDED_BY,
    NODE_PUBLICATIONS,
    NODE_SYNONYMS,
    NODE_TAXON,
    NODE_UNMATCHED_NAME,
    NODE_URLS,
    PROJECT_ROOT,
)

DEFAULT_RANKING_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "prefix_ranking.yaml"

_UNRANKED = 10_000

# Babel types an id by its structure, and many drugs have no Drug-typed id: tezepelumab is ChemicalEntity from Babel,
# and the aggregators' Drug + Protein spans two families, so its other ids inherit nothing. ATC gives a 7-character
# code almost only to a therapeutic substance, so a chemical holding one is a drug whatever its ids are typed. Set here,
# after clustering, it cannot change what merges. Chemicals only: the same codes sit on organisms (probiotics), cells,
# devices and proteins. A few 7-character codes name a group instead; these are the ones found on chemical nodes.
DRUG_CATEGORY = "biolink:Drug"
ATC_SUBSTANCE_CODE = re.compile(r"ATC:[A-Z]\d{2}[A-Z]{2}\d{2}")
ATC_GROUP_CODES = frozenset(
    {
        "ATC:A07XA03",  # calcium compounds
        "ATC:B05BA03",  # carbohydrates (on the ChEBI class "carbohydrate")
        "ATC:D11AC02",  # cadmium compounds
        "ATC:D11AC08",  # sulfur compounds
        "ATC:S01AX03",  # zinc compounds
    }
)


@dataclass(frozen=True)
class NameCorrection:
    """A chemical's display name, checked by hand against its structure. It applies to the cluster holding the audited
    structure and the id it is keyed by, while that id's own label is still the right name."""

    wrong_name: str
    name: str
    inchikey: str


# Id whose own label is the right name -> correction. Each wrong name is an HMDB or aggregator label for another
# compound (#29). Keyed by that id, not by the representative, which a new member can change. The wrong name leaves
# the synonyms too, so it is matched exactly: the same string is a fine name on the compound it belongs to.
GLYCOPEPTIDE = "(alpha-D-mannosyl)7-beta-D-mannosyl-diacetylchitobiosyl-L-asparagine, isoform {} (protein)"
NAME_CORRECTIONS: dict[str, NameCorrection] = {
    "CHEBI:37675": NameCorrection(
        wrong_name="Glycoprotein-phospho-D-mannose", name="aldehydo-D-mannose", inchikey="GZCGUPFRVQAUEE-KVTDHHQDSA-N"
    ),
    "CHEBI:27931": NameCorrection(
        wrong_name="Amylose", name="alpha-maltotriose", inchikey="FYGDTMLNYKFZSV-PXXRMHSHSA-N"
    ),
    "CHEBI:16449": NameCorrection(
        wrong_name=GLYCOPEPTIDE.format("A"), name="alanine", inchikey="QNAYBMKLOCPYGJ-UHFFFAOYSA-N"
    ),
    "MESH:C108415": NameCorrection(
        wrong_name=GLYCOPEPTIDE.format("B"), name="4-iodoanisole", inchikey="SYSZENVIJHPFNL-UHFFFAOYSA-N"
    ),
}


def _name_correction(cluster_ids: set[str], label_of: Callable[[str], str | None]) -> NameCorrection | None:
    for curie, correction in NAME_CORRECTIONS.items():
        if (
            curie in cluster_ids
            and f"INCHIKEY:{correction.inchikey}" in cluster_ids
            and label_of(curie) == correction.name
        ):
            return correction
    return None


class PrefixRanking:
    """Ranked prefixes per family, for representative selection."""

    def __init__(self, rankings: dict[str, list[str]], default: list[str], last: Sequence[str] = ()):
        self.rankings = rankings
        self.default = default
        self.last = {prefix: _UNRANKED + 1 + i for i, prefix in enumerate(last)}

    @classmethod
    def load(cls, path: str | Path | None = None) -> PrefixRanking:
        path = Path(path) if path is not None else DEFAULT_RANKING_PATH
        data = yaml.safe_load(Path(path).read_text()) or {}
        return cls(data.get("rankings", {}), data.get("_default", []), data.get("_last", []))

    def order_for(self, family: str | None) -> dict[str, int]:
        prefixes = self.rankings.get(family) if family else None
        if not prefixes:
            prefixes = self.default
        return {prefix: i for i, prefix in enumerate(prefixes)} | self.last

    def rank_key(self, curie: str, order: Mapping[str, int]) -> tuple[int, str]:
        prefix = curie.split(":", 1)[0]
        return (order.get(prefix, _UNRANKED), curie)


def _cluster_family(members: Sequence[dict], families: BranchFamilies) -> str | None:
    """Pick the family whose ranking should govern the representative.

    Uses the branch intersection across members (the guardrail-valid family set);
    on an all-wildcard cluster returns None (default ranking).
    """
    node_branches = [families.branches(m.get(NODE_CATEGORIES)) for m in members]
    branches = BranchFamilies.cluster_branches(node_branches)
    if branches is None:
        # a still-conflicting cluster (repair left it intact): fall back to the
        # union so we at least pick a sensible ranking
        union: set[str] = set()
        for nb in node_branches:
            if nb is not ALL_FAMILIES and nb != ALL_FAMILIES:
                union |= nb
        branches = frozenset(union)
    if branches is ALL_FAMILIES or branches == ALL_FAMILIES or not branches:
        return None
    return sorted(branches)[0]


def _union_field(members: Iterable[dict], key: str) -> list[str]:
    out: set[str] = set()
    for m in members:
        val = m.get(key)
        if isinstance(val, list):
            out.update(v for v in val if v)
        elif val:
            out.add(val)
    return sorted(out)


def _highest_ranked_matching(
    members: Sequence[dict],
    order: Mapping[str, int],
    ranking: PrefixRanking,
    predicate: Callable[[dict], bool],
) -> dict | None:
    """The member of highest rank satisfying ``predicate`` (deterministic)."""
    candidates = [m for m in members if predicate(m)]
    if not candidates:
        return None
    return min(candidates, key=lambda m: ranking.rank_key(m[NODE_ID], order))


def materialize_cluster(
    members: Sequence[dict],
    ranking: PrefixRanking,
    families: BranchFamilies,
    label_of: Callable[[str], str | None] = lambda _curie: None,
) -> dict:
    """Reconcile member node dicts into one canonical node (order-independent).

    ``label_of`` returns Babel's per-id label for a CURIE (or None); the
    display name prefers the REPRESENTATIVE id's own Babel label (attributable to that
    exact canonical id), falling back to the highest-ranked member name when Babel has
    none, and to the highest-ranked unmatched name when no member has a name. Every other
    member name is retained as a synonym either way."""
    if not members:
        raise ValueError("cannot materialize an empty cluster")

    family = _cluster_family(members, families)
    order = ranking.order_for(family)

    representative = min(members, key=lambda m: ranking.rank_key(m[NODE_ID], order))
    name_src = _highest_ranked_matching(members, order, ranking, lambda m: bool(m.get(NODE_NAME)))
    unmatched_name_src = _highest_ranked_matching(members, order, ranking, lambda m: bool(m.get(NODE_UNMATCHED_NAME)))
    desc_src = _highest_ranked_matching(members, order, ranking, lambda m: bool(m.get(NODE_DESCRIPTION)))

    # union of all member ids and their equivalent id lists
    equivalent_ids = set(_union_field(members, NODE_EQUIVALENT_IDS))
    equivalent_ids.update(m[NODE_ID] for m in members)

    # synonyms = union of member synonym BAGS, but only from members whose own
    # equivalency list stays WITHIN this cluster. A member that "spills" (its equiv
    # list names ids the clustering split into OTHER clusters -- e.g. an aggregator
    # node that fused 1300 Reactome ids into a gene) is conflated relative to us, and
    # its bag is not per-id attributable, so we keep only its primary name (below),
    # not the bag. Members with no/subset equiv lists are self-consistent -> trusted.
    cluster_ids = {m[NODE_ID] for m in members}
    synonyms: set[str] = set()
    for m in members:
        member_equiv = {e for e in (m.get(NODE_EQUIVALENT_IDS) or ()) if e}
        if member_equiv <= cluster_ids:
            synonyms.update(s for s in (m.get(NODE_SYNONYMS) or ()) if s)
    # Display name: the representative id's OWN Babel label wins (per-id, attributable);
    # else the highest-ranked member name. Every other member name (incl. the rep's
    # own source name when the Babel label differs) is retained as a synonym below.
    # Where no id has a name, a name its source kept out of the matching: an ICD code that is to merge into its
    # disease, when the disease has no room for it (ICD10CM:I21.9, where I21 already sat in acute myocardial
    # infarction, had no name at all).
    chosen_name = (
        label_of(representative[NODE_ID])
        or (name_src[NODE_NAME] if name_src else None)
        or (unmatched_name_src[NODE_UNMATCHED_NAME] if unmatched_name_src else None)
    )
    correction = _name_correction(cluster_ids, label_of)
    if correction and chosen_name == correction.wrong_name:
        logging.info(
            "materialize: %s is named %r, not %r (hand-checked)",
            representative[NODE_ID],
            correction.name,
            correction.wrong_name,
        )
        chosen_name = correction.name
    for m in members:
        # keep each member's source name AND its own Babel label (only an audited wrong name goes, below)
        for candidate in (m.get(NODE_NAME), m.get(NODE_UNMATCHED_NAME), label_of(m[NODE_ID])):
            if candidate and candidate != chosen_name:
                synonyms.add(candidate)
    synonyms.discard(chosen_name)
    if correction:
        synonyms.discard(correction.wrong_name)

    # categories: leaf-filtered union is deferred to the caller (needs BiolinkClient);
    # here we emit the raw union (plus a derived Drug) and let resolve/postprocess leaf-filter.
    categories = _union_field(members, NODE_CATEGORIES)
    if family == CHEMICAL_FAMILY and any(
        ATC_SUBSTANCE_CODE.fullmatch(curie) and curie not in ATC_GROUP_CODES for curie in cluster_ids
    ):
        categories = sorted({*categories, DRUG_CATEGORY})

    # taxon: guardrail guarantees <=1 distinct; None on any residual conflict
    taxa = {m.get(NODE_TAXON) for m in members if m.get(NODE_TAXON)}
    taxon = next(iter(taxa)) if len(taxa) == 1 else None

    canonical: dict = {
        NODE_ID: representative[NODE_ID],
        NODE_CATEGORIES: categories,
        NODE_PROVIDED_BY: _union_field(members, NODE_PROVIDED_BY),
        NODE_EQUIVALENT_IDS: sorted(equivalent_ids),
    }
    if chosen_name:
        canonical[NODE_NAME] = chosen_name
    if desc_src and desc_src.get(NODE_DESCRIPTION):
        canonical[NODE_DESCRIPTION] = desc_src[NODE_DESCRIPTION]
    if synonyms:
        canonical[NODE_SYNONYMS] = sorted(synonyms)
    if taxon:
        canonical[NODE_TAXON] = taxon

    urls = _union_field(members, NODE_URLS)
    if urls:
        canonical[NODE_URLS] = urls
    publications = _union_field(members, NODE_PUBLICATIONS)
    if publications:
        canonical[NODE_PUBLICATIONS] = publications

    # scalar chemistry props: take from the highest-ranked member that has them
    for key in (NODE_CHEMICAL_FORMULA, NODE_EXACT_MASS):
        src = _highest_ranked_matching(members, order, ranking, lambda m, k=key: m.get(k) not in (None, ""))
        if src is not None:
            canonical[key] = src[key]

    # attributes: union of second-level keys (deterministic; keep earliest by rank)
    attributes = _merge_attributes(members, order, ranking)
    if attributes:
        canonical[NODE_ATTRIBUTES] = attributes

    return canonical


def _merge_attributes(members: Sequence[dict], order: Mapping[str, int], ranking: PrefixRanking) -> dict:
    merged: dict = {}
    for m in sorted(members, key=lambda m: ranking.rank_key(m[NODE_ID], order)):
        attrs = m.get(NODE_ATTRIBUTES)
        if not isinstance(attrs, dict):
            continue
        for k, v in attrs.items():
            merged.setdefault(k, v)  # highest-ranked member wins per key
    return merged
