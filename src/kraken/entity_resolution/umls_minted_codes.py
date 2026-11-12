"""Put an id UMLS minted itself into the node of the concept it was minted in, when nothing else ties it down.

OMIM has no code for its clinical synopsis terms, so UMLS mints one: OMIM:MTHU069260 is the term "Type 2 diabetes"
of concept UMLS:C0011860. The code means nothing outside that concept: kg2 carries 9,060 of them, each with exactly
one concept. But kg2 reports that as a UMLS:xref close_match, below tau, and the term's name often differs from the
concept's ("Type 2 diabetes" / "Diabetes Mellitus, Non-Insulin-Dependent"), so 6,591 of them shipped in 2.4.0 as
one-id nodes. One of those carries real edges: 10 Mendelian randomization edges from the multiomics KGs that the
type 2 diabetes node never got.

This runs AFTER clustering, rather than as heavier evidence, so it cannot move anything else. Measured over the
44k ids these codes can reach, when clustering was label propagation, heavier evidence tore 115 clusters and fused
469; even limited to the codes with no other tie it tore 87, because a new leaf changed how the rest clustered. A
code that already shares a node stays there: its other tie disagrees with the concept, which is not a call to make
here.

A code may only add itself to the node. OMIM ranks above HP, so in 2,262 of the 2.4.0 nodes the code would have
become the node's id, and in 2,258 it would have added a type the node's own ids do not have (kg2 types it from the
concept's UMLS semantic type). And where the node's id has no Babel label, the code's name could become the node's:
318 nodes would have been renamed. Those codes stay apart. The ones that join lend the node no description either
(see ``_stage4_materialize``).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Callable, Mapping

from kraken.entity_resolution.families import ALL_FAMILIES, CHEMICAL_FAMILY, BranchFamilies
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo, cluster_violations
from kraken.entity_resolution.materialize import PrefixRanking

UMLS_XREF_RELATION = "UMLS:xref"
UMLS_MINTED_PREFIX = "OMIM:MTH"


def is_umls_minted_code(curie: str, relation: str) -> bool:
    """Whether kg2's original edge ``curie -> relation -> ...`` says which UMLS concept a minted code belongs to."""
    return relation == UMLS_XREF_RELATION and curie.startswith(UMLS_MINTED_PREFIX)


def attach_orphan_codes(
    code_to_concept: Mapping[str, str],
    curie_to_cluster: dict[str, int],
    *,
    info_of: Callable[[str], NodeInfo],
    categories_of: Mapping[str, tuple[str, ...]],
    label_of: Callable[[str], str | None],
    ranking: PrefixRanking,
    guardrail_config: GuardrailConfig,
    next_cluster_id: int,
) -> tuple[set[str], int]:
    """Move each code that is alone into its concept's cluster, in place, if the guardrails allow it and the node
    keeps its id, name, types and taxon. Returns the codes moved and the next free cluster id."""
    wanted = {curie_to_cluster[c] for pair in code_to_concept.items() for c in pair if c in curie_to_cluster}
    members_of: dict[int, list[str]] = defaultdict(list)
    for curie, cluster in curie_to_cluster.items():
        if cluster in wanted:
            members_of[cluster].append(curie)

    attached: set[str] = set()
    refused = changes_node = 0
    for code, concept in sorted(code_to_concept.items()):
        own_cluster = curie_to_cluster.get(code)
        if own_cluster is not None and len(members_of[own_cluster]) > 1:
            continue
        cluster = curie_to_cluster.get(concept)
        node = members_of[cluster] if cluster is not None else [concept]
        ids = sorted({*node, code})
        info = {curie: info_of(curie) for curie in ids}
        # An untyped code would be a guardrail wildcard, so nothing could refuse it.
        if info[code].branches is ALL_FAMILIES or cluster_violations(ids, info, guardrail_config):
            refused += 1
            continue
        # The family orders the members, and so decides the node's id and which member describes it. The name stays
        # only if it is that id's Babel label: without one, materialize falls back to the best-ranked member's name,
        # which the code could be. A chemical's name can be another spelling its members vote on, so none join one.
        # A clinical term has no species: a code with a taxon is not one of these, and its records may disagree on it
        # beyond the one taxon per id the guardrails see, which materialize would then drop from the node.
        family = _family(node, info)
        leader = _leader(node, family, ranking)
        node_categories = {category for curie in node for category in categories_of.get(curie, ())}
        keeps_node = (
            family != CHEMICAL_FAMILY
            and _family(ids, info) == family
            and _leader(ids, family, ranking) == leader
            and label_of(leader)
            and set(categories_of.get(code, ())) <= node_categories
            and info[code].taxon is None
        )
        if not keeps_node:
            changes_node += 1
            continue
        if cluster is None:
            cluster, next_cluster_id = next_cluster_id, next_cluster_id + 1
            curie_to_cluster[concept] = cluster
            members_of[cluster].append(concept)
        curie_to_cluster[code] = cluster
        members_of[cluster].append(code)
        attached.add(code)
    logging.info(
        "entity_resolution: put %d UMLS-minted codes into their concept's node; %d refused by the guardrails, %d "
        "would have changed the node's id, name, types or taxon",
        len(attached),
        refused,
        changes_node,
    )
    return attached, next_cluster_id


def _family(ids: list[str], info: Mapping[str, NodeInfo]) -> str | None:
    """The family whose ranking materialize orders the node's members by."""
    branches = BranchFamilies.cluster_branches([info[curie].branches for curie in ids])
    return None if not branches or branches == ALL_FAMILIES else sorted(branches)[0]


def _leader(ids: list[str], family: str | None, ranking: PrefixRanking) -> str:
    """The id materialize will make the node's id: the top-ranked one for its family."""
    order = ranking.order_for(family)
    return min(ids, key=lambda curie: ranking.rank_key(curie, order))
