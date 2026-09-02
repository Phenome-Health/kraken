"""Clustering: connected components, then Leiden/CPM per non-trivial component.

Pipeline (plan §2):

1. **Connected components** (union-find) — decomposition + parallelism; trivial
   components (isolated pairs handled, singletons added by ``resolve``) need no
   Leiden.
2. **Leiden with CPM** on each non-trivial component — the primary mechanism. It
   is what separates ``Adams-Oliver syndrome 1`` from ``AOS2``, which no
   guardrail can see. CPM avoids modularity's resolution limit.

``gamma``'s plain reading: A and B joined by total weight w stay separate iff
``gamma * |A| * |B| > w``.

**Determinism is required** (releases are DOI-archived): fixed seed and sorted
node order everywhere.
"""

from __future__ import annotations

import logging
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
        # path compression
        while self.parent[x] != root:
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


def leiden_cpm(
    nodes: list[str],
    edges: list[WeightedPair],
    gamma: float,
    *,
    seed: int = DEFAULT_SEED,
) -> list[list[str]]:
    """Run Leiden with CPM on one component. Deterministic given seed + sorted
    node order. Returns sub-clusters as sorted node lists.

    ``nodes`` must contain every endpoint in ``edges``.
    """
    import igraph  # imported lazily so the module loads without the C deps
    import leidenalg

    ordered = sorted(nodes)
    index = {name: i for i, name in enumerate(ordered)}
    g = igraph.Graph()
    g.add_vertices(len(ordered))
    ig_edges = [(index[a], index[b]) for a, b, _w in edges]
    ig_weights = [w for _a, _b, w in edges]
    g.add_edges(ig_edges)
    if ig_weights:
        g.es["weight"] = ig_weights
    partition = leidenalg.find_partition(
        g,
        leidenalg.CPMVertexPartition,
        weights="weight" if ig_weights else None,
        resolution_parameter=gamma,
        n_iterations=-1,
        seed=seed,
    )
    clusters = [sorted(ordered[i] for i in community) for community in partition]
    clusters.sort(key=lambda members: members[0])
    return clusters


def cluster_pairs(
    pairs: Iterable[WeightedPair],
    gamma: float,
    *,
    seed: int = DEFAULT_SEED,
) -> list[list[str]]:
    """Full clustering of the match graph: components then Leiden/CPM.

    Returns every multi-node and single-node cluster arising from ``pairs``.
    (CURIEs never appearing in ``pairs`` are added as singletons by ``resolve``.)
    """
    pair_list = list(pairs)
    components = connected_components(pair_list)

    # index edges by component for the non-trivial ones
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
        elif len(comp) == 2:
            # CPM on an isolated pair degenerates to: split iff gamma > w.
            (_a, _b, w) = comp_edges[ci][0]
            if w > gamma:
                clusters.append(comp)
            else:
                clusters.append([comp[0]])
                clusters.append([comp[1]])
        else:
            clusters.extend(leiden_cpm(comp, comp_edges[ci], gamma, seed=seed))
    logging.info(
        "clustering: %d components -> %d clusters from %d pairs",
        len(components),
        len(clusters),
        len(pair_list),
    )
    return clusters
