"""Put a WHO ICD-10 code into the node of the ICD-10-CM code it is.

The icd10 source says which ICD-10-CM code a WHO code is: the one with the same code and title, 5,575 of them
(see ``harmonizers/icd10.py``).

This runs AFTER clustering, rather than as match evidence, so that it moves the WHO code and nothing else. As
evidence the statement weighs 1, more than the name match that holds an ICD-10-CM code in most nodes, so the two
merged first. They then stood as a group of two ids, which needs more than tau per id to join a group of several,
and one name match of 0.45 is not that: in the 2.4.0-dev5 build 1,015 ICD-10-CM codes lost every id they had
shared a node with (``ICD10CM:B18.0`` left NCIT's, SNOMED's and UMLS's concept of its name).

The code is its twin, so the guardrails see it as an id of its twin's type, and it lends the node no type,
description or name of its own (see ``_stage4_materialize``): its names are synonyms there, and name the node only
when no other id has a name of any kind. Nor the id: ``ICD10`` ranks after every other prefix. Nor a species: a
code some source gives one stays where it is, as an ICD code has none.

kg2 and ROBOKOP list 68 of these codes on a node, 55 of them on another node than their ICD-10-CM code's and seldom
the right one (``ICD10:S12.8`` "Fracture of other parts of neck" on Mondo's tracheal disorder). Such a code leaves
that node for its twin's, and the list's statement stays in the graph as a close_match edge between the two nodes.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping

from kraken.entity_resolution.final_clusters import FinalClusters
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo, cluster_violations
from kraken.utils.constants import SAME_AS_PREDICATE

ICD10_SOURCE = "icd10"


def is_twin_statement(source: str, predicate: str) -> bool:
    """Whether a harmonized edge of ``source`` says its subject, a WHO code, is the ICD-10-CM code that is its
    object. The icd10 source states no other equivalence."""
    return source == ICD10_SOURCE and predicate == SAME_AS_PREDICATE


def attach_twins(
    twin_of: Mapping[str, str],
    curie_to_cluster: dict[str, int],
    *,
    info_of: Callable[[str], NodeInfo],
    guardrail_config: GuardrailConfig,
    next_cluster_id: int,
) -> tuple[set[str], int]:
    """Move each WHO code that has no taxon into the cluster of its ICD-10-CM code, in place, if the guardrails
    allow it there as an id of that code's type. Returns the codes moved and the next free cluster id.

    Which codes move does not depend on their order. A code is checked, and counts from then on, as an id of its
    twin's type with no taxon. Its twin is in the node, so that breaks no rule on families or taxa, whichever other
    WHO codes are in the node, come or about to leave. (A one-id rule on ``ICD10`` itself could refuse it; none is
    set.)"""
    clusters = FinalClusters(curie_to_cluster, {*twin_of, *twin_of.values()}, next_cluster_id)
    moved: dict[str, NodeInfo] = {}  # each with what it counts as from then on
    there_already = from_another_node = refused = 0
    for code, twin in sorted(twin_of.items()):
        node = clusters.members_with(twin)
        if code in node:
            there_already += 1
            continue
        info = {curie: moved.get(curie) or info_of(curie) for curie in node}
        info[code] = NodeInfo(curie=code, branches=info[twin].branches)
        if info_of(code).taxon is not None or cluster_violations([*node, code], info, guardrail_config):
            refused += 1
            continue
        from_another_node += not clusters.is_alone(code)
        clusters.join(code, twin)
        moved[code] = info[code]
    logging.info(
        "entity_resolution: put %d WHO ICD-10 codes into the node of their ICD-10-CM code, %d of them out of a node "
        "another source had them on; %d were there already, %d have a taxon or were refused by the guardrails",
        len(moved),
        from_another_node,
        there_already,
        refused,
    )
    return set(moved), clusters.next_cluster_id
