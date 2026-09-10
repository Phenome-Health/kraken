"""Clustering: connected components, then label propagation per component.

Pipeline:

1. **Connected components** (union-find) — decomposition + parallelism; trivial
   components need no work.
2. **Label propagation** on each non-trivial component (weighted). Unlike
   Leiden/CPM there is no resolution parameter: label propagation merges what is
   connected, and the guardrails (branch / taxon / one-id) do the splitting. The
   match graph is already shaped for this — evidence below tau is dropped, and
   guardrail-conflicting edges are pruned before clustering — so LP labels the
   surviving connected structure and the guardrails break up anything that
   crossed a hard boundary transitively.

**Determinism** (releases are DOI-archived): fixed seed (igraph RNG seeded) and
sorted node order everywhere.
"""

from __future__ import annotations

import logging
import random
from collections.abc import Iterable

from kraken.entity_resolution.match_graph import WeightedPair

DEFAULT_SEED = 20240101


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def add(self, x: str) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x: str) -> str:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:  # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def connected_components(pairs: Iterable[WeightedPair]) -> list[list[str]]:
    """Partition the pair graph into connected components.

    Returns components as sorted node lists, ordered deterministically by their
    smallest member. Only nodes appearing in ``pairs`` are included.
    """
    uf = _UnionFind()
    for a, b, _w in pairs:
        uf.add(a)
        uf.add(b)
        uf.union(a, b)
    comps: dict[str, list[str]] = {}
    for node in uf.parent:
        comps.setdefault(uf.find(node), []).append(node)
    result = [sorted(members) for members in comps.values()]
    result.sort(key=lambda members: members[0])
    return result


def label_propagation(
    nodes: list[str],
    edges: list[WeightedPair],
    *,
    seed: int = DEFAULT_SEED,
) -> list[list[str]]:
    """Weighted label propagation on one component. Deterministic given seed +
    sorted node order. Returns sub-clusters as sorted node lists.

    ``nodes`` must contain every endpoint in ``edges``. With no edges, every node
    is its own cluster.
    """
    import igraph  # imported lazily so the module loads without the C deps

    ordered = sorted(nodes)
    index = {name: i for i, name in enumerate(ordered)}
    g = igraph.Graph()
    g.add_vertices(len(ordered))
    g.add_edges([(index[a], index[b]) for a, b, _w in edges])
    weights = [w for _a, _b, w in edges]

    # Seed igraph's RNG so label propagation is reproducible run-to-run.
    try:
        igraph.set_random_number_generator(random.Random(seed))
    except Exception:  # pragma: no cover - older/newer igraph API differences
        random.seed(seed)
    partition = g.community_label_propagation(weights=weights or None)

    clusters = [sorted(ordered[i] for i in community) for community in partition]
    clusters.sort(key=lambda members: members[0])
    return clusters


def cluster_pairs(
    pairs: Iterable[WeightedPair],
    *,
    seed: int = DEFAULT_SEED,
) -> list[list[str]]:
    """Full clustering of the match graph: components then label propagation.

    Returns every multi-node and single-node cluster arising from ``pairs``.
    (CURIEs never appearing in ``pairs`` are added as singletons by the caller.)
    """
    pair_list = list(pairs)
    components = connected_components(pair_list)

    node_to_comp: dict[str, int] = {}
    for ci, comp in enumerate(components):
        for node in comp:
            node_to_comp[node] = ci
    comp_edges: dict[int, list[WeightedPair]] = {i: [] for i in range(len(components))}
    for a, b, w in pair_list:
        comp_edges[node_to_comp[a]].append((a, b, w))

    clusters: list[list[str]] = []
    for ci, comp in enumerate(components):
        if len(comp) == 1:
            clusters.append(comp)
        else:
            clusters.extend(label_propagation(comp, comp_edges[ci], seed=seed))
    logging.info(
        "clustering: %d components -> %d clusters from %d pairs",
        len(components),
        len(clusters),
        len(pair_list),
    )
    return clusters
