"""Materialize a cluster into one canonical node (plan §1, §7 representative).

Order-independence is required (the whole point of the overhaul), so every
reconciliation here is a function of the *set* of members, not their order:

* union fields (equivalent_ids, synonyms, provided_by, urls, publications) are
  set-unions;
* the representative id / name / description are chosen by a **total order**
  (ranked prefix for the cluster's family, then CURIE), which is
  order-independent;
* categories are the leaf-filtered union;
* taxon is the single guardrail-permitted taxon (or None on conflict).

The representative id is cosmetic (queries resolve through the membership map),
so it only labels the cluster; name and description may come from different
members (a lower-ranked id can supply the description the top id lacks).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path

import yaml

from kraken.entity_resolution.families import ALL_FAMILIES, BranchFamilies
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
    NODE_URLS,
    PROJECT_ROOT,
)

DEFAULT_RANKING_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "prefix_ranking.yaml"

_UNRANKED = 10_000


class PrefixRanking:
    """Ranked prefixes per family, for representative selection."""

    def __init__(self, rankings: dict[str, list[str]], default: list[str]):
        self.rankings = rankings
        self.default = default

    @classmethod
    def load(cls, path: str | Path | None = None) -> PrefixRanking:
        path = Path(path) if path is not None else DEFAULT_RANKING_PATH
        data = yaml.safe_load(Path(path).read_text()) or {}
        return cls(data.get("rankings", {}), data.get("_default", []))

    def order_for(self, family: str | None) -> dict[str, int]:
        prefixes = self.rankings.get(family) if family else None
        if not prefixes:
            prefixes = self.default
        return {prefix: i for i, prefix in enumerate(prefixes)}

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
) -> dict:
    """Reconcile member node dicts into one canonical node (order-independent)."""
    if not members:
        raise ValueError("cannot materialize an empty cluster")

    family = _cluster_family(members, families)
    order = ranking.order_for(family)

    representative = min(members, key=lambda m: ranking.rank_key(m[NODE_ID], order))
    name_src = _highest_ranked_matching(members, order, ranking, lambda m: bool(m.get(NODE_NAME)))
    desc_src = _highest_ranked_matching(members, order, ranking, lambda m: bool(m.get(NODE_DESCRIPTION)))

    # union of all member ids and their equivalent id lists
    equivalent_ids = set(_union_field(members, NODE_EQUIVALENT_IDS))
    equivalent_ids.update(m[NODE_ID] for m in members)

    # synonyms = union of member synonyms plus any member names not chosen as THE name
    synonyms = set(_union_field(members, NODE_SYNONYMS))
    chosen_name = name_src.get(NODE_NAME) if name_src else None
    for m in members:
        nm = m.get(NODE_NAME)
        if nm and nm != chosen_name:
            synonyms.add(nm)
    synonyms.discard(chosen_name)

    # categories: leaf-filtered union is deferred to the caller (needs BiolinkClient);
    # here we emit the raw union and let resolve/postprocess leaf-filter.
    categories = _union_field(members, NODE_CATEGORIES)

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
