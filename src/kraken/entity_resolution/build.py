"""Memory-bounded entity resolution over the full harmonized graph.

Produces the two artifacts the build needs:

* the integrated **canonical nodes** file (one node per cluster), and
* a ``node_id -> representative_curie`` map for resolving edge endpoints,

replacing ``integrate.py``'s legacy equivalency-trusting node merge.

The build must stay under ~48 GB RAM at ~30M CURIEs, so the two structures that
would otherwise dominate memory are kept **out of core**:

* **evidence accumulation** — all evidence is streamed to a temp file and
  combined with an external ``sort`` (the pattern ``integrate_edges`` already
  uses), so peak memory is one CURIE-pair's evidence, not the whole graph;
* **node materialization** — harmonized nodes are streamed, tagged with their
  cluster, external-sorted by cluster, and reconciled one cluster at a time.

The numeric core (factorization, connected components) uses numpy/scipy, which
are both far faster and far more memory-compact than Python dicts/objects:
CURIE strings become int codes, components come from ``scipy`` csgraph, and
Leiden runs per non-trivial component so peak clustering memory is bounded by the
largest component rather than the whole graph.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import jsonlines
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from kraken.entity_resolution.clustering import DEFAULT_SEED, leiden_cpm
from kraken.entity_resolution.families import BranchFamilies
from kraken.entity_resolution.guardrails import (
    GuardrailConfig,
    NodeInfo,
    cluster_violations,
    enforce_cluster,
    ids_per_cluster_histogram,
    log_oversized_clusters,
)
from kraken.entity_resolution.match_graph import clique_evidence, match_predicate_evidence
from kraken.entity_resolution.materialize import PrefixRanking, materialize_cluster
from kraken.entity_resolution.name_sim import DEFAULT_STOPLIST, is_droppable, normalize_name
from kraken.entity_resolution.sri_nodenorm import NodeNormClient
from kraken.entity_resolution.weights import NAME_SIMILARITY_GROUP, ERWeights
from kraken.utils.constants import (
    EDGE_ATTRIBUTES,
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_SUBJECT,
    KG2_INFORES,
    NODE_CATEGORIES,
    NODE_EQUIVALENT_IDS,
    NODE_ID,
    NODE_NAME,
    NODE_PROVIDED_BY,
    NODE_TAXON,
)
from kraken.utils.kg_io import remove_file, stream_edges_from_jsonl, stream_nodes_from_jsonl

SEP = "\t"

# Aggregators that pre-merge via SRI NN / Babel: their harmonized edge subject/object
# are canonicalized ids, so a match-predicate edge must be un-canonicalized to its
# ORIGINAL endpoints before entering the match graph (else it just re-imports Babel).
# KG2 keeps the originals in ``kg2pre_ids`` (handled below). The other four emit no
# match-predicate edges today (verified), so we have no format to parse for them yet;
# their match edges are skipped until they do / their format is known.
CANONICALIZED_AGGREGATOR_SOURCES: frozenset[str] = frozenset(
    {"kg2", "robokop", "translator-kg-open", "microbiome-kg", "multiomics-kg"}
)
# KG2's per-edge original ids: "orig_subject---relation---q---q---q---orig_object---src".
KG2_PRE_IDS_ATTR = "kg2pre_ids"
_KG2_ID_SEP = "---"


def _original_endpoints(edge: dict, source: str) -> list[tuple[str, str]] | None:
    """Original (pre-canonicalization) subject/object pairs for a match-predicate edge.

    A Babel-canonicalizing aggregator stores merged endpoints; the originals live in edge
    attributes, and one merged edge can carry several original pairs. KG2: parse them out
    of ``kg2pre_ids`` (subject at field 0, object at field 5). Native (non-canonicalized)
    sources just use their own subject/object. Returns None for a canonicalized aggregator
    we can't yet un-canonicalize, so its match edges are skipped rather than ingested as
    Babel's own clustering.
    """
    if source == "kg2":
        attrs = (edge.get(EDGE_ATTRIBUTES) or {}).get(KG2_INFORES) or {}
        pairs: list[tuple[str, str]] = []
        for raw in attrs.get(KG2_PRE_IDS_ATTR) or []:
            parts = raw.split(_KG2_ID_SEP)
            if len(parts) >= 6 and parts[0] not in ("", "None") and parts[5] not in ("", "None"):
                pairs.append((parts[0], parts[5]))
        return pairs or None
    if source in CANONICALIZED_AGGREGATOR_SOURCES:
        return None
    return [(edge.get(EDGE_SUBJECT, ""), edge.get(EDGE_OBJECT, ""))]


def _external_sort(input_path: Path, output_path: Path, key_args: list[str], temp_dir: Path) -> None:
    """Byte-ordered (deterministic) external sort with spill files on ``temp_dir``."""
    subprocess.run(
        ["sort", "-t", SEP, *key_args, "-T", str(temp_dir), "-o", str(output_path), str(input_path)],
        check=True,
        env={**os.environ, "LC_ALL": "C"},
    )


# --------------------------------------------------------------------------------------
# Stage 1: stream harmonized data -> evidence file, name file, per-CURIE guardrail facts
# --------------------------------------------------------------------------------------


def _stage1_write_evidence_and_facts(
    config,
    weights: ERWeights,
    evidence_path: Path,
    names_path: Path,
) -> tuple[dict[str, set[str]], dict[str, str], set[str]]:
    """One streaming pass over harmonized nodes+edges. Writes equivalency-clique and
    match-predicate evidence, writes ``normalized_name<TAB>curie`` rows for name
    similarity, and returns per-CURIE categories/taxon plus the set of harmonized
    node primary ids (so bare equivalency-list ids can be told apart later).

    Category strings are interned so the facts dict stays compact (there are only
    ~150 distinct Biolink categories).
    """
    node_cats: dict[str, set[str]] = {}
    node_taxon: dict[str, str] = {}
    node_ids: set[str] = set()

    with open(evidence_path, "w") as ev, open(names_path, "w") as nm:
        for source, (nodes_path, edges_path) in sorted(config.all_harmonized_paths_resolved.items()):
            if Path(nodes_path).exists():
                for node in stream_nodes_from_jsonl(Path(nodes_path)):
                    node_id = node.get(NODE_ID)
                    if not node_id:
                        continue
                    node_ids.add(node_id)
                    for a, b, group, weight in clique_evidence(node.get(NODE_EQUIVALENT_IDS) or [], source, weights):
                        ev.write(f"{a}{SEP}{b}{SEP}{group}{SEP}{weight}\n")
                    cats = node.get(NODE_CATEGORIES) or []
                    if cats:
                        node_cats.setdefault(node_id, set()).update(sys.intern(c) for c in cats)
                    taxon = node.get(NODE_TAXON)
                    if taxon:
                        node_taxon[node_id] = sys.intern(taxon)
                    name_norm = normalize_name(node.get(NODE_NAME))
                    if not is_droppable(name_norm, min_length=weights.min_name_length, stoplist=DEFAULT_STOPLIST):
                        nm.write(f"{name_norm}{SEP}{node_id}\n")
            # Match-predicate (close/exact/same_as) edges are match-graph evidence, but
            # only on their ORIGINAL endpoints (see _original_endpoints): a Babel-
            # canonicalized endpoint would just re-import Babel's clustering. KG2 also
            # needs its close_match down-weighted where subclass edges co-occur, so it
            # gets a dedicated two-phase writer.
            if Path(edges_path).exists():
                if source == "kg2":
                    _write_kg2_match_evidence(Path(edges_path), weights, ev)
                else:
                    for edge in stream_edges_from_jsonl(Path(edges_path)):
                        predicate = edge.get(EDGE_PREDICATE, "")
                        if weights.predicate_weight(predicate) is None:
                            continue  # not a usable match predicate; skip before un-canon work
                        endpoints = _original_endpoints(edge, source)
                        if endpoints is None:
                            continue  # canonicalized aggregator without a known un-canonicalizer
                        for subject, object_ in endpoints:
                            ev_edge = match_predicate_evidence(subject, object_, predicate, source, weights)
                            if ev_edge is not None:
                                a, b, group, weight = ev_edge
                                ev.write(f"{a}{SEP}{b}{SEP}{group}{SEP}{weight}\n")
    return node_cats, node_taxon, node_ids


# subclass_of / superclass_of between the same pair signals the co-occurring close_match
# is a mislabeled hierarchical relation, not equivalence, so we down-weight it.
SUBCLASS_PREDICATES: frozenset[str] = frozenset({"biolink:subclass_of", "biolink:superclass_of"})


def _subclass_penalized_weight(base_weight: float, hierarchical_count: int, decay: float) -> float:
    """Down-weight a close_match by ``decay`` for each co-occurring hierarchical edge
    (more hierarchical evidence -> weaker close_match)."""
    return base_weight * (decay**hierarchical_count)


def _write_kg2_match_evidence(edges_path: Path, weights: ERWeights, ev) -> None:
    """Emit KG2 match-predicate evidence on un-canonicalized endpoints, down-weighting
    each close_match by how many subclass/superclass edges the same original pair has.

    Two-phase over KG2's edges (option (a)): first count hierarchical edges per original
    pair and buffer the match pairs (bounded by KG2's edge count, not the whole graph),
    then emit each with its penalty applied.
    """
    hierarchical_counts: dict[tuple[str, str], int] = defaultdict(int)
    match_pairs: list[tuple[str, str, str]] = []  # (a, b, predicate), a <= b
    for edge in stream_edges_from_jsonl(edges_path):
        predicate = edge.get(EDGE_PREDICATE, "")
        is_hierarchical = predicate in SUBCLASS_PREDICATES
        if not is_hierarchical and weights.predicate_weight(predicate) is None:
            continue
        for subject, object_ in _original_endpoints(edge, "kg2") or []:
            if subject == object_:
                continue
            a, b = (subject, object_) if subject <= object_ else (object_, subject)
            if is_hierarchical:
                hierarchical_counts[(a, b)] += 1
            else:
                match_pairs.append((a, b, predicate))

    group = weights.correlation_group("kg2")
    for a, b, predicate in match_pairs:
        base = weights.predicate_weight(predicate)
        weight = _subclass_penalized_weight(base, hierarchical_counts.get((a, b), 0), weights.subclass_penalty_decay)
        ev.write(f"{a}{SEP}{b}{SEP}{group}{SEP}{weight}\n")


def _stage1b_append_name_similarity(names_path: Path, evidence_path: Path, weights: ERWeights, temp_dir: Path) -> None:
    """Group CURIEs by normalized name (external sort) and append name-similarity
    clique evidence for each group within the size cap. Bounded memory: one name
    group at a time."""
    sorted_names = temp_dir / "er_s1_names_sorted.tmp"
    try:
        _external_sort(names_path, sorted_names, ["-k1,1"], temp_dir)
        with open(sorted_names) as fin, open(evidence_path, "a") as ev:
            current = None
            ids: list[str] = []

            def flush(group_ids: list[str]) -> None:
                unique_ids = sorted(set(group_ids))
                if not (2 <= len(unique_ids) <= weights.name_group_cap):
                    return
                w = weights.name_similarity_weight
                for i in range(len(unique_ids)):
                    for j in range(i + 1, len(unique_ids)):
                        a, b = unique_ids[i], unique_ids[j]
                        if a > b:
                            a, b = b, a
                        ev.write(f"{a}{SEP}{b}{SEP}{NAME_SIMILARITY_GROUP}{SEP}{w}\n")

            for line in fin:
                name, _, curie = line.rstrip("\n").partition(SEP)
                if name != current and ids:
                    flush(ids)
                    ids = []
                current = name
                ids.append(curie)
            if ids:
                flush(ids)
    finally:
        remove_file(sorted_names)


# --------------------------------------------------------------------------------------
# Stage 2: accumulate evidence -> tau-filtered weighted pairs (external sort)
# --------------------------------------------------------------------------------------


def _stage2_accumulate_pairs(evidence_path: Path, pairs_path: Path, weights: ERWeights, temp_dir: Path) -> int:
    """Combine evidence per CURIE pair (max within correlation group, sum across)
    and keep pairs meeting tau. Returns the number of pairs written."""
    sorted_ev = temp_dir / "er_s2_evidence_sorted.tmp"
    n_pairs = 0
    try:
        _external_sort(evidence_path, sorted_ev, ["-k1,1", "-k2,2", "-k3,3"], temp_dir)
        with open(sorted_ev) as fin, open(pairs_path, "w") as out:
            cur_a = cur_b = None
            group_max: dict[str, float] = {}

            def flush() -> int:
                if cur_a is None:
                    return 0
                total = sum(group_max.values())
                if total >= weights.tau:
                    out.write(f"{cur_a}{SEP}{cur_b}{SEP}{total}\n")
                    return 1
                return 0

            for line in fin:
                a, b, group, weight_s = line.rstrip("\n").split(SEP)
                weight = float(weight_s)
                if a != cur_a or b != cur_b:
                    n_pairs += flush()
                    cur_a, cur_b = a, b
                    group_max = {}
                prev = group_max.get(group)
                if prev is None or weight > prev:
                    group_max[group] = weight
            n_pairs += flush()
    finally:
        remove_file(sorted_ev)
    return n_pairs


# --------------------------------------------------------------------------------------
# Stage 3: components + Leiden + guardrails -> curie -> cluster id
# --------------------------------------------------------------------------------------


def _stage3_cluster(
    pairs_path: Path,
    weights: ERWeights,
    families: BranchFamilies,
    guardrail_config: GuardrailConfig,
    harmonized_cats: dict[str, set[str]],
    node_taxon: dict[str, str],
    node_ids: set[str],
    nodenorm: NodeNormClient,
    seed: int,
) -> tuple[dict[str, int], dict[str, dict[int, int]], dict[str, str], dict[str, tuple[str, ...]]]:
    """Cluster the weighted pair graph. Returns ``curie -> cluster_id`` for every
    CURIE appearing in a pair, the ids-per-prefix histogram, ``curie -> label`` for
    bare ids, and ``curie -> single intrinsic category`` for every match-graph node.

    Uses int codes (numpy) and scipy connected components; Leiden + guardrails run
    per non-trivial component so peak memory is bounded by the largest component.
    """
    if os.path.getsize(pairs_path) == 0:
        return {}, {}, {}, {}

    df = pd.read_csv(pairs_path, sep=SEP, names=["a", "b", "w"], dtype={"a": str, "b": str, "w": "float32"})
    codes, uniques = pd.factorize(pd.concat([df["a"], df["b"]], ignore_index=True), sort=False)
    n = len(df)
    a_codes = codes[:n].astype(np.int64)
    b_codes = codes[n:].astype(np.int64)
    w = df["w"].to_numpy()
    num_nodes = len(uniques)
    del df

    # Connected components on the undirected graph.
    graph = coo_matrix((w, (a_codes, b_codes)), shape=(num_nodes, num_nodes))
    n_components, labels = connected_components(graph, directed=False)
    logging.info("entity_resolution: %d nodes, %d pairs, %d connected components", num_nodes, n, n_components)

    # Edges grouped by component (both endpoints share a label): sort edge indices
    # by component so each component's edges are a contiguous slice.
    edge_labels = labels[a_codes]
    edge_order = np.argsort(edge_labels, kind="stable")
    edge_labels_sorted = edge_labels[edge_order]
    comp_edge_starts = np.searchsorted(edge_labels_sorted, np.arange(n_components), side="left")
    comp_edge_ends = np.searchsorted(edge_labels_sorted, np.arange(n_components), side="right")

    # Node codes grouped by component.
    node_order = np.argsort(labels, kind="stable")
    labels_sorted = labels[node_order]
    comp_node_starts = np.searchsorted(labels_sorted, np.arange(n_components), side="left")
    comp_node_ends = np.searchsorted(labels_sorted, np.arange(n_components), side="right")

    # ONE intrinsic category per match-graph node. A match-graph node is a single
    # unmerged identifier, so it must carry exactly one category (its own), NOT the
    # possibly-conflated category list of a harmonized node it appears in. The node
    # normalizer is the source of truth (it types each id individually); on a miss
    # we fall back to a lone harmonized category if the source gave exactly one,
    # else NamedThing (empty -> wildcard). This is what lets the branch guardrail
    # bite: a conflated MONDO node becomes just Disease, not {Disease,Gene,Protein}.
    all_curies = [uniques[i] for i in range(num_nodes)]
    logging.info("entity_resolution: resolving %d match-graph node categories via node normalizer", len(all_curies))
    resolved = nodenorm.resolve(all_curies)
    mg_categories: dict[str, tuple[str, ...]] = {}
    bare_names: dict[str, str] = {}
    for curie in all_curies:
        norm = resolved.get(curie)
        cats = tuple(norm.categories) if norm and norm.categories else ()
        if not cats:
            harm = harmonized_cats.get(curie)
            if harm and len(harm) == 1:  # trust a lone source category; ignore conflated multi
                cats = tuple(harm)
        mg_categories[curie] = cats
        if norm and norm.label and curie not in node_ids:
            bare_names[curie] = norm.label

    def info_provider(curie: str) -> NodeInfo:
        return NodeInfo(curie=curie, categories=mg_categories.get(curie, ()), taxon=node_taxon.get(curie))

    curie_to_cluster: dict[str, int] = {}
    all_clusters: list[list[str]] = []
    next_cluster_id = 0

    for comp in range(n_components):
        node_idx = node_order[comp_node_starts[comp] : comp_node_ends[comp]]
        member_curies = [uniques[c] for c in node_idx]
        if len(member_curies) == 1:
            raw_clusters = [member_curies]
        else:
            edge_idx = edge_order[comp_edge_starts[comp] : comp_edge_ends[comp]]
            comp_edges = [(uniques[a_codes[i]], uniques[b_codes[i]], float(w[i])) for i in edge_idx]
            info = {c: info_provider(c) for c in member_curies}
            # Guardrail-aware FORMATION: drop edges between conflicting nodes before
            # clustering, so incompatible things never merge in the first place (this
            # can also split the component for free). cluster_violations on the pair
            # covers all enforced guardrails: branch, one_id, taxon.
            comp_edges = [
                e for e in comp_edges if not cluster_violations([e[0], e[1]], info, families, guardrail_config)
            ]
            if len(member_curies) == 2:
                merged = bool(comp_edges) and comp_edges[0][2] > weights.gamma
                raw_clusters = [member_curies] if merged else [[member_curies[0]], [member_curies[1]]]
            else:
                raw_clusters = leiden_cpm(member_curies, comp_edges, weights.gamma, seed=seed)
            # Guardrails as the backstop, split until valid (catches transitive
            # conflicts the pairwise prune can't see).
            adjacency = _adjacency(comp_edges)
            splitter = _make_splitter(adjacency, weights.gamma, seed)
            checked: list[list[str]] = []
            for cluster in raw_clusters:
                checked.extend(
                    enforce_cluster(cluster, info, families, guardrail_config, splitter=splitter, adjacency=adjacency)
                )
            raw_clusters = checked

        for cluster in raw_clusters:
            for curie in cluster:
                curie_to_cluster[curie] = next_cluster_id
            all_clusters.append(cluster)
            next_cluster_id += 1

    log_oversized_clusters(all_clusters, guardrail_config)
    return curie_to_cluster, ids_per_cluster_histogram(all_clusters), bare_names, mg_categories


def _adjacency(edges: list[tuple[str, str, float]]) -> dict[str, dict[str, float]]:
    adj: dict[str, dict[str, float]] = defaultdict(dict)
    for a, b, weight in edges:
        adj[a][b] = weight
        adj[b][a] = weight
    return adj


def _make_splitter(adjacency, gamma, seed):
    def splitter(members: list[str]) -> list[list[str]]:
        member_set = set(members)
        seen: set[tuple[str, str]] = set()
        edges: list[tuple[str, str, float]] = []
        for a in members:
            for b, weight in adjacency.get(a, {}).items():
                if b in member_set:
                    key = (a, b) if a < b else (b, a)
                    if key not in seen:
                        seen.add(key)
                        edges.append((key[0], key[1], weight))
        if not edges:
            return [members]
        return leiden_cpm(members, edges, gamma * 2.0, seed=seed)

    return splitter


# --------------------------------------------------------------------------------------
# Stage 4: materialize canonical nodes (external group-by) -> node_id -> representative
# --------------------------------------------------------------------------------------


def _stage4_materialize(
    config,
    curie_to_cluster: dict[str, int],
    mg_categories: dict[str, tuple[str, ...]],
    bare_names: dict[str, str],
    node_ids: set[str],
    ranking: PrefixRanking,
    families: BranchFamilies,
    biolink,
    temp_dir: Path,
) -> dict[str, str]:
    """Stream harmonized nodes, group by cluster on disk, reconcile one cluster at
    a time, and write the canonical nodes file. Returns ``node_id -> representative``.

    Materialization is cluster-driven: bare member ids (from canonicalized sources'
    equivalency lists, with no harmonized node of their own) are emitted as
    synthetic member dicts carrying their node-normalizer name/category, so an
    entity that splits off an aggregator clique becomes its own named node instead
    of being dropped. A harmonized node whose id is not in any pair is its own
    singleton cluster (keyed by its id). Peak memory is one cluster's members.
    """
    keyed = temp_dir / "er_s4_nodes_keyed.tmp"
    keyed_sorted = temp_dir / "er_s4_nodes_keyed_sorted.tmp"
    node_id_to_rep: dict[str, str] = {}

    # Inverse of curie_to_cluster: the canonical node's equivalent_ids must be its
    # exact cluster membership, NOT the union of member source lists (which can
    # contain ids the guardrails split into other clusters, breaking disjointness).
    cluster_members: dict[int, list[str]] = defaultdict(list)
    for curie, cid in curie_to_cluster.items():
        cluster_members[cid].append(curie)

    try:
        with jsonlines.open(keyed, "w") as writer:
            for _source, (nodes_path, _edges) in sorted(config.all_harmonized_paths_resolved.items()):
                if not Path(nodes_path).exists():
                    continue
                for node in stream_nodes_from_jsonl(Path(nodes_path)):
                    node_id = node.get(NODE_ID)
                    if not node_id:
                        continue
                    # Replace the (possibly conflated) source category list with this
                    # id's single intrinsic category, so the merged node's categories
                    # are the union of its members' true types, not conflation leftovers.
                    node[NODE_CATEGORIES] = sorted(mg_categories.get(node_id, ()))
                    cid = curie_to_cluster.get(node_id)
                    key = f"c{cid}" if cid is not None else f"s:{node_id}"
                    writer.write([key, node])
            # Synthetic member dicts for bare ids (no harmonized node of their own).
            for curie, cid in curie_to_cluster.items():
                if curie in node_ids:
                    continue
                synthetic = {
                    NODE_ID: curie,
                    NODE_CATEGORIES: sorted(mg_categories.get(curie, ())),
                    NODE_PROVIDED_BY: [],
                }
                label = bare_names.get(curie)
                if label:
                    synthetic[NODE_NAME] = label
                writer.write([f"c{cid}", synthetic])
        _external_sort(keyed, keyed_sorted, ["-k1,1"], temp_dir)

        with (
            jsonlines.open(keyed_sorted) as reader,
            jsonlines.open(config.integrated_nodes_path, "w") as out,
        ):
            current_key: str | None = None
            members: list[dict] = []

            def flush(group_key: str, group: list[dict]) -> None:
                if not group:
                    return
                node = materialize_cluster(group, ranking, families)
                if biolink is not None and node.get(NODE_CATEGORIES):
                    leaves = biolink.filter_to_leaf_categories(node[NODE_CATEGORIES])
                    if leaves:
                        node[NODE_CATEGORIES] = sorted(leaves)
                # equivalent_ids = exact cluster membership (guarantees disjointness)
                if group_key.startswith("c"):
                    node[NODE_EQUIVALENT_IDS] = sorted(cluster_members[int(group_key[1:])])
                else:  # singleton: the node's own id(s)
                    node[NODE_EQUIVALENT_IDS] = sorted({m[NODE_ID] for m in group})
                out.write(node)
                rep = node[NODE_ID]
                for member in group:
                    node_id_to_rep[member[NODE_ID]] = rep

            for key, node in reader:
                if key != current_key and members:
                    flush(current_key, members)
                    members = []
                current_key = key
                members.append(node)
            if current_key is not None:
                flush(current_key, members)
    finally:
        remove_file(keyed)
        remove_file(keyed_sorted)
    return node_id_to_rep


# --------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------


def resolve_entities(config, biolink) -> dict[str, str]:
    """Run entity resolution end to end (out of core) and write the canonical nodes
    file. Returns ``node_id -> representative_curie`` for edge resolution."""
    weights = ERWeights.load()
    families = BranchFamilies.load()
    ranking = PrefixRanking.load()
    guardrail_config = GuardrailConfig()
    temp_dir = config.integrated_dir
    temp_dir.mkdir(parents=True, exist_ok=True)

    evidence_path = temp_dir / "er_s1_evidence.tmp"
    names_path = temp_dir / "er_s1_names.tmp"
    pairs_path = temp_dir / "er_s2_pairs.tmp"

    config.er_nodenorm_cache_path.parent.mkdir(parents=True, exist_ok=True)
    nodenorm = NodeNormClient(config.er_nodenorm_cache_path)
    try:
        logging.info("entity_resolution: streaming evidence + facts...")
        harmonized_cats, node_taxon, node_ids = _stage1_write_evidence_and_facts(
            config, weights, evidence_path, names_path
        )
        _stage1b_append_name_similarity(names_path, evidence_path, weights, temp_dir)

        logging.info("entity_resolution: accumulating pairs...")
        _stage2_accumulate_pairs(evidence_path, pairs_path, weights, temp_dir)

        logging.info("entity_resolution: clustering...")
        curie_to_cluster, histogram, bare_names, mg_categories = _stage3_cluster(
            pairs_path,
            weights,
            families,
            guardrail_config,
            harmonized_cats,
            node_taxon,
            node_ids,
            nodenorm,
            DEFAULT_SEED,
        )
        _log_histogram(histogram)

        logging.info("entity_resolution: materializing canonical nodes...")
        node_id_to_rep = _stage4_materialize(
            config, curie_to_cluster, mg_categories, bare_names, node_ids, ranking, families, biolink, temp_dir
        )
    finally:
        nodenorm.close()
        remove_file(evidence_path)
        remove_file(names_path)
        remove_file(pairs_path)

    _report_eval(curie_to_cluster)
    logging.info("entity_resolution: %d node ids mapped to representatives", len(node_id_to_rep))
    return node_id_to_rep


def _log_histogram(histogram: dict[str, dict[int, int]]) -> None:
    for prefix in ("HGNC", "NCBIGene", "MONDO", "UniProtKB", "RM", "LM"):
        counts = histogram.get(prefix)
        if counts:
            multi = {k: v for k, v in sorted(counts.items()) if k > 1}
            if multi:
                logging.info("entity_resolution ids-per-cluster[%s]: >1 -> %s", prefix, multi)


def _report_eval(curie_to_cluster: dict[str, int]) -> None:
    from kraken.entity_resolution.eval.scorer import DEFAULT_GROUND_TRUTH_PATH, load_gold, score

    if not DEFAULT_GROUND_TRUTH_PATH.exists() or not curie_to_cluster:
        return
    res = score(load_gold(DEFAULT_GROUND_TRUTH_PATH), curie_to_cluster)
    logging.info(
        "entity_resolution eval: precision=%.4f recall=%.4f f1=%.4f (ml covered=%d, cl covered=%d)",
        res.precision,
        res.recall,
        res.f1,
        res.must_link_covered,
        res.cannot_link_covered,
    )
