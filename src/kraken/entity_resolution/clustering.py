"""Clustering: one connected component of the match graph in, its clusters out.

The clusters are a function of the evidence only -- the pairs and their weights, the guardrail facts, the Babel
cliques -- and never of id spelling, input order or float rounding: releases are DOI-archived, and the benchmark has
to be able to measure a small change. So weights are integers, strengths are compared exactly, and links of equal
strength are decided together:

1. Each Babel clique the guardrails allow starts as one group; every other id starts alone.
2. The link between two groups is the summed weight of the pairs between them. Its strength is that weight per id
   of the smaller group.
3. The strongest links are decided first, all links of one strength at once, on the groups as they stood before any
   of them was applied:

   * a link between groups the guardrails keep apart is dropped;
   * a link merges its groups when it has tau per id of the smaller one (see ``strong_enough``), and a weaker one
     waits: a later merge can add weight to it;
   * a group with equally strong links to groups that cannot all be one group is in an exact tie: its links to
     the ones that rival each other are not applied;
   * links that chain several groups are applied as one union, or not at all when the guardrails forbid the union.

   The links of a tie and of a refused chain are withdrawn for good: equal claims that contradict each other
   cancel, and only the remaining evidence can still merge those groups. Ties that held to the end and refused
   chains are reported (see ``Tie``).
"""

from __future__ import annotations

import heapq
import json
import logging
from collections import Counter, defaultdict
from collections.abc import Collection, Iterable, Mapping, Sequence
from fractions import Fraction
from functools import reduce
from pathlib import Path
from typing import NamedTuple

from kraken.entity_resolution.guardrails import GroupState, GuardrailConfig, NodeInfoMap, group_state, join_states
from kraken.entity_resolution.weights import WEIGHT_UNIT

SEP = "\t"


class Tie(NamedTuple):
    """Equally strong links that merged nothing, as reported for curation.

    ``tie``: ``group`` has such links to the ``candidates``, which cannot all be one group with it.
    ``chain``: such links chain ``group`` and the ``candidates``, and the guardrails forbid them as one group.
    Two that end as the same row are reported once: a row is about the groups as they ended.
    """

    kind: str
    group: str  # groups are named by their smallest id, and sized, as the component ended
    size: int
    strength: float  # of each of those links: weight per id of the smaller group
    candidates: tuple[str, ...]


def cluster_component(
    members: Sequence[str],
    pairs: Sequence[tuple[str, str, int]],
    info: NodeInfoMap,
    config: GuardrailConfig,
    clique_of: Mapping[str, int],
    tau: int,
) -> tuple[list[list[str]], list[Tie]]:
    """The clusters of one component, and its ties and refused chains (see the module docstring).

    ``members`` are a whole component of the above-tau pairs, formed before any pair was pruned: it holds every
    starting member of each of its cliques, so a clique is judged as a whole. ``pairs`` are its pairs that no
    guardrail forbids, and ``clique_of`` gives the Babel clique of each id that starts in one. Weights and ``tau``
    are in ``WEIGHT_UNIT``.
    """
    if len(members) == 2 and len(pairs) == 1:  # most components: one pair, nothing to decide
        return [sorted(members)], []

    # A group is known by the position of one of its members, and no decision reads that number. ``name``, a
    # group's smallest id, only labels the report.
    parent = list(range(len(members)))
    size = [1] * len(members)
    name = list(members)
    state = [group_state(info[member], config.enforced_prefixes) for member in members]
    link: list[dict[int, int]] = [{} for _ in members]  # link[a][b] == link[b][a]: the weight between two groups

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def add_link(a: int, b: int, weight: int) -> None:
        link[a][b] = link[b][a] = link[a].get(b, 0) + weight

    def drop_link(a: int, b: int) -> None:
        del link[a][b], link[b][a]

    def merge(groups: Collection[int]) -> list[tuple[int, int]] | None:
        """Make the groups one, or return None when the guardrails forbid the union. Returns the links whose
        strength the merge changed."""
        joined = reduce(join_states, (state[group] for group in groups))
        if joined is None:
            return None
        keep = max(groups, key=size.__getitem__)
        size_before = size[keep]
        took_weight = set()
        for gone in groups:
            if gone == keep:
                continue
            parent[gone] = keep
            size[keep] += size[gone]
            name[keep] = min(name[keep], name[gone])
            for other, weight in link[gone].items():
                del link[other][gone]
                if other != keep:
                    add_link(keep, other, weight)
                    took_weight.add(other)
            link[gone] = {}
        state[keep] = joined
        # A strength is weight per id of the smaller group: the other links of ``keep`` still have theirs. Deciding
        # every link again after every merge is what would make a dense component quadratic.
        return [(keep, other) for other in link[keep] if other in took_weight or size[other] > size_before]

    in_clique: dict[int, list[int]] = defaultdict(list)
    for i, member in enumerate(members):
        if member in clique_of:
            in_clique[clique_of[member]].append(i)
    for clique in in_clique.values():
        merge(clique)  # a clique that breaks a guardrail by itself starts as single ids
    group_of = {member: find(i) for i, member in enumerate(members)}
    for a, b, weight in pairs:
        if group_of[a] != group_of[b]:
            add_link(group_of[a], group_of[b], weight)

    def strength(a: int, b: int) -> Fraction:
        return Fraction(link[a][b], min(size[a], size[b]))  # exact: no float decides which links are equal

    def strong_enough(weight: int, smaller: int) -> bool:
        """Whether a link of this weight merges two groups, the smaller of this size. Reaching tau is enough for a
        single id, as it is for a pair. Two groups of several ids must exceed it: at exactly tau, two name matches
        of 0.45 merged insulin with the three ids of a proinsulin drug."""
        return weight >= tau * smaller if smaller == 1 else weight > tau * smaller

    def rivalled(group: int, partners: list[int]) -> list[int]:
        """The partners (the groups ``group`` has equally strong links to) that have a rival: other partners such
        that they cannot all be one group with ``group``."""
        # Guardrails. They read only a group's state, so they are asked once per kind of partner, not once per pair
        # of partners: a hub has hundreds of partners of one kind.
        kinds = Counter(state[partner] for partner in partners)

        def kept_apart(kind: GroupState) -> bool:
            with_group = join_states(state[group], kind)
            return any(join_states(with_group, other) is None for other in kinds if other != kind or kinds[kind] > 1)

        apart = {kind: kept_apart(kind) for kind in kinds}
        # Evidence, two partners at a time. The three are one group only if, once ``group`` has joined one partner,
        # the other still merges: its link to ``group`` plus the pairs between the two partners must be
        # strong_enough for the smaller side. A partner whose own link is strong_enough for its whole size can
        # always come second. The rest are heavy, and two heavy partners are rivals unless the pairs between them
        # make up for it.
        heavy = {partner for partner in partners if not strong_enough(link[group][partner], size[partner])}

        def has_heavy_rival(partner: int) -> bool:
            made_up = sum(
                strong_enough(link[group][other] + between, min(size[partner], size[other]))
                for other, between in link[partner].items()
                if other in heavy
            )
            return made_up < len(heavy) - 1

        rivals = [p for p in partners if apart[state[p]] or (p in heavy and has_heavy_rival(p))]
        # Taxa and enforced prefixes clash two at a time, but families can clash three at a time. So when the
        # partners left still cannot all join ``group``, they share no family, and each typed one is a rival.
        left = set(partners).difference(rivals)
        if reduce(join_states, {state[partner] for partner in left}, state[group]) is None:
            rivals += [partner for partner in left if state[partner][0] is not None]
        return rivals

    # heapq pops the smallest entry, so the strengths are negated.
    heap = [(-strength(a, b), a, b) for a, links in enumerate(link) for b in links if a < b]
    heapq.heapify(heap)
    withdrawn: list[tuple[str, Fraction, list[int]]] = []  # kind, strength, groups (a tie's own group first)
    while heap:
        negated = heap[0][0]
        level = -negated
        current: set[tuple[int, int]] = set()  # the links that have this strength now: heap entries go stale
        while heap and heap[0][0] == negated:
            _negated, a, b = heapq.heappop(heap)
            if b in link[a] and strength(a, b) == level:
                current.add((a, b))
        mergeable = []
        for a, b in current:
            if join_states(state[a], state[b]) is None:
                drop_link(a, b)  # for good: what a guardrail forbids two groups, it forbids all they grow into
            elif strong_enough(link[a][b], min(size[a], size[b])):
                mergeable.append((a, b))
        partners: dict[int, list[int]] = defaultdict(list)
        for a, b in mergeable:
            partners[a].append(b)
            partners[b].append(a)
        tied: set[tuple[int, int]] = set()
        for group, around in partners.items():
            rivals = rivalled(group, around) if len(around) > 1 else []
            if rivals:
                tied.update((min(group, rival), max(group, rival)) for rival in rivals)
                withdrawn.append(("tie", level, [group, *rivals]))
        for a, b in tied:
            drop_link(a, b)
        changed = []
        for chain in _chains([pair for pair in mergeable if pair not in tied]):
            groups = {group for pair in chain for group in pair}
            links = merge(groups)
            if links is None:
                for a, b in chain:
                    drop_link(a, b)
                withdrawn.append(("chain", level, sorted(groups)))
            else:
                changed += links
        for a, b in changed:
            if b in link[a]:  # a later chain of this level may have merged ``b`` away, and listed the link itself
                heapq.heappush(heap, (-strength(a, b), min(a, b), max(a, b)))

    clusters: dict[int, list[str]] = defaultdict(list)
    for i, member in enumerate(members):
        clusters[find(i)].append(member)
    reports = set()
    for kind, level, groups in withdrawn:
        own, *others = (find(group) for group in groups)
        if kind == "chain":
            own, *others = sorted({own, *others}, key=name.__getitem__)
        elif own in others:
            continue  # the tie did not hold: other links took the group to one of its candidates after all
        candidates = tuple(sorted({name[other] for other in others}))
        reports.add(Tie(kind, name[own], size[own], float(level) / WEIGHT_UNIT, candidates))
    return sorted(sorted(cluster) for cluster in clusters.values()), sorted(reports)


def _chains(links: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """The links grouped into chains: the links of one chain connect its groups, directly or through each other."""
    top: dict[int, int] = {}

    def find(group: int) -> int:
        while top.setdefault(group, group) != group:
            top[group] = top[top[group]]
            group = top[group]
        return group

    for a, b in links:
        top[find(a)] = find(b)
    chains: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for pair in links:
        chains[find(pair[0])].append(pair)
    return list(chains.values())


def log_ties(ties: Iterable[Tie], report_path: Path) -> None:
    """Write the ties and refused chains for curation: one usually means two sources disagree, or a guardrail fact
    is wrong. The candidates are a JSON list, because an id can hold a space (``HCPCS:Level 1: A0000-A9999``)."""
    rows = sorted(ties)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w") as out:
        out.write(f"kind{SEP}group{SEP}size{SEP}strength{SEP}candidates\n")
        for tie in rows:
            candidates = json.dumps(tie.candidates)
            out.write(f"{tie.kind}{SEP}{tie.group}{SEP}{tie.size}{SEP}{tie.strength}{SEP}{candidates}\n")
    chains = sum(tie.kind == "chain" for tie in rows)
    logging.info(
        "entity_resolution: %d exact ties and %d refused chains merged nothing; listed in %s",
        len(rows) - chains,
        chains,
        report_path,
    )
