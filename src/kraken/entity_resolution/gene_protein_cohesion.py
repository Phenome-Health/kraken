"""Put back the clusters a single Babel gene/protein clique was split across.

A Babel clique within ``clique_cap`` is a FULL clique -- every member vouching for every other at the Babel weight
-- and a gene's clique and its protein's clique are conflated into one (see ``build._write_babel_evidence``). Label
propagation is a community algorithm, though, not a connectivity one, and a conflated clique is two dense halves
joined at that one uniform weight. When one half has outside evidence and the other has none, the halves settle
into different labels and the clique tears in two.

ACE is the shape of it. All 11 ids sit in one Babel clique, so each has 10 edges of 0.5, and no guardrail objects.
But the 4 gene ids also carry 6-9 points of NCBI Gene, OMIM, ENSEMBL and LOINC evidence, so they hold the gene
label against 3.5 of Babel pull, while the 7 protein ids (PR, UniProtKB, MESH, NCIT, UMLS) have no outside evidence
at all and only compare Babel against Babel: 3.0 among themselves against 2.0 to the gene side. The bigger half
wins, and the gene and its protein end up in different nodes -- which our own policy says is wrong, since a gene
and the protein it encodes in one species are one entity.

So: when a Babel clique's members land in several clusters and EVERY id of those clusters is a gene or protein (or
untyped), the clusters are rejoined -- if the guardrails allow the union. The guardrails are what keep this honest:
measured on the 2.3.0 build, about 43,500 gene/protein cliques span more than one cluster, and only about a quarter
of them are rejoined here. The other three quarters are refused on TAXON, and rightly: they are Babel cliques that
mix orthologs across species (rat RGD ids sitting in a dog's or a pig's clique), which must stay apart.

The repair is deliberately limited to this one family. Tried across all families it is a wash -- on the benchmark
it trades 420 must-link pairs recovered against 397 cannot-link pairs lost -- because elsewhere a torn Babel clique
is usually label propagation doing real work, quietly separating a class from its members and a part from its
whole (it would put HMDB's "Glycolipids" back with sphingomyelin, undoing ``babel_outliers``, and merge
mesonephros with adult mammalian kidney, neurocranium with chondrocranium, a GO process with a LOINC part).
Restricted to gene/protein it costs nothing measurable: no cannot-link pair in the benchmark is lost.

Babel clique outliers need no special handling here: ``babel_outliers`` never judges gene/protein ids, whose names
are symbols rather than descriptions.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Callable, Iterator
from pathlib import Path

from kraken.entity_resolution.families import ALL_FAMILIES
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo, cluster_violations

SEP = "\t"
GENE_PROTEIN_FAMILY = "gene_protein"

InfoProvider = Callable[[str], NodeInfo]


def _cliques(cliques_path: Path) -> Iterator[list[str]]:
    """Each clique's members, from ``id<TAB>hub<TAB>size`` rows written one clique at a time."""
    members: list[str] = []
    current = None
    with open(cliques_path) as fin:
        for line in fin:
            curie, hub, _size = line.rstrip("\n").split(SEP)
            if hub != current:
                if members:
                    yield members
                current, members = hub, []
            members.append(curie)
    if members:
        yield members


def _gene_protein(info: NodeInfo) -> bool:
    """A gene or protein -- or untyped, which is a guardrail wildcard and can be carried along."""
    return info.branches is ALL_FAMILIES or set(info.branches) == {GENE_PROTEIN_FAMILY}


def rejoin_split_gene_protein_cliques(
    curie_to_cluster: dict[str, int],
    cliques_path: Path,
    info_of: InfoProvider,
    guardrail_config: GuardrailConfig,
) -> list[list[list[str]]]:
    """Rejoin the clusters that one Babel gene/protein clique was split across, in place.

    Returns the clusters merged, as one list of member lists per merge, so the caller can retally what it counted
    per cluster. A group whose union breaks a guardrail is left alone entirely -- taxon is the usual reason, and a
    clique that mixes two species has no valid subset worth guessing at.
    """
    # 1. The clusters one clique's gene/protein members were split across, unioned as we stream.
    parent: dict[int, int] = {}

    def find(cluster: int) -> int:
        while parent[cluster] != cluster:
            parent[cluster] = parent[parent[cluster]]
            cluster = parent[cluster]
        return cluster

    def union(one: int, other: int) -> None:
        for cluster in (one, other):
            parent.setdefault(cluster, cluster)
        root_one, root_other = find(one), find(other)
        if root_one != root_other:
            parent[max(root_one, root_other)] = min(root_one, root_other)

    split_cliques = 0
    for members in _cliques(cliques_path):
        clustered = [m for m in members if m in curie_to_cluster]
        clusters = {curie_to_cluster[m] for m in clustered}
        if len(clusters) < 2 or not all(_gene_protein(info_of(m)) for m in clustered):
            continue
        split_cliques += 1
        first, *rest = sorted(clusters)
        for other in rest:
            union(first, other)

    if not parent:
        logging.info("entity_resolution: no Babel gene/protein clique was split across clusters")
        return []

    # 2. Those clusters' FULL membership -- the guardrails judge the whole union, not just the clique.
    roots = {cluster: find(cluster) for cluster in parent}
    members_of: dict[int, list[str]] = defaultdict(list)
    for curie, cluster in curie_to_cluster.items():
        if cluster in roots:
            members_of[cluster].append(curie)

    # 3. Merge each group the guardrails allow.
    by_root: dict[int, list[int]] = defaultdict(list)
    for cluster, root in roots.items():
        by_root[root].append(cluster)
    merged: list[list[list[str]]] = []
    refused = 0
    for root, clusters in by_root.items():
        ids = [curie for cluster in clusters for curie in members_of[cluster]]
        info = {curie: info_of(curie) for curie in ids}
        if not all(_gene_protein(info[curie]) for curie in ids) or cluster_violations(
            sorted(ids), info, guardrail_config
        ):
            refused += 1
            continue
        for curie in ids:
            curie_to_cluster[curie] = root
        merged.append([members_of[cluster] for cluster in sorted(clusters)])

    logging.info(
        "entity_resolution: %d Babel gene/protein cliques were split across clusters; rejoined %d groups "
        "(%d ids), left %d alone because the union breaks a guardrail (usually two species in one clique)",
        split_cliques,
        len(merged),
        sum(len(cluster) for group in merged for cluster in group),
        refused,
    )
    return merged
