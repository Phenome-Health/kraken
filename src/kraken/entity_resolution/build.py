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
CURIE strings become int codes, components come from ``scipy`` csgraph, and label
propagation runs per non-trivial component so peak clustering memory is bounded by
the largest component rather than the whole graph.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import jsonlines
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from kraken.entity_resolution.clustering import DEFAULT_SEED, label_propagation
from kraken.entity_resolution.families import ALL_FAMILIES, BranchFamilies
from kraken.entity_resolution.guardrails import (
    GuardrailConfig,
    NodeInfo,
    cluster_violations,
    enforce_cluster,
    ids_per_cluster_histogram,
    log_oversized_clusters,
)
from kraken.entity_resolution.match_graph import alias_evidence, clique_evidence, match_predicate_evidence
from kraken.entity_resolution.materialize import PrefixRanking, materialize_cluster
from kraken.entity_resolution.name_sim import DEFAULT_STOPLIST, is_droppable, normalize_name
from kraken.entity_resolution.sri_nodenorm import NodeNormClient, infer_category, infer_taxon
from kraken.entity_resolution.uncanonicalize import (
    ALIAS_EVIDENCE_SOURCES,
    CANONICALIZED_AGGREGATOR_SOURCES,
    is_coarser_than_canonical,
    original_alias_pairs,
)
from kraken.entity_resolution.uncanonicalize import (
    original_endpoints as _original_endpoints,
)
from kraken.entity_resolution.weights import NAME_SIMILARITY_GROUP, ERWeights
from kraken.utils.constants import (
    EDGE_PREDICATE,
    EDGE_PRIMARY_KS,
    NODE_CATEGORIES,
    NODE_EQUIVALENT_IDS,
    NODE_ID,
    NODE_NAME,
    NODE_PROVIDED_BY,
    NODE_TAXON,
)
from kraken.utils.kg_io import remove_file, stream_edges_from_jsonl, stream_nodes_from_jsonl

SEP = "\t"

# Set in an id's ``seeds`` mask when it was seeded ONLY by an alias (not by any node or equivalency list), so the
# edge pass can tell "already in the graph" from "introduced by an earlier alias". Far above any source's bit.
ALIAS_SEED_BIT = 1 << 62


def _evidence_row(a: str, b: str, group: str, weight: float, *, kind: str) -> str:
    """One evidence line. ``kind`` (e.g. "equiv:kg2", "alias:robokop", "nn", "name_sim") is not used to weigh
    anything -- ``group`` does that -- but survives to the oversized-cluster diagnostics, which is the only way to
    tell a normalizer clique from an aggregator list once both are in the sri_nn_derived group."""
    return f"{a}{SEP}{b}{SEP}{group}{SEP}{weight}{SEP}{kind}\n"


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
    families: BranchFamilies,
    evidence_path: Path,
    names_path: Path,
    source_bits: dict[str, int],
) -> tuple[dict[str, set[str]], dict[str, str], set[str], dict[str, int], dict[str, set[str]]]:
    """One streaming pass over harmonized nodes+edges. Writes equivalency-clique
    (native sources only) and match-predicate evidence, writes
    ``normalized_name<TAB>curie`` rows for name similarity, and returns:

    * ``inherited_cats``: curie -> candidate categories, propagated from every node
      whose categories are cleanly SINGLE-family onto each id in that node's
      equivalency list. This is how the vast majority of ids (which only ever appear
      as equiv-list members, never as a primary node) get a category. Conflated
      multi-family nodes are deliberately skipped so their mis-typing doesn't spread.
    * per-CURIE harmonized taxon, the set of harmonized node primary ids, and
      ``seeds`` (every curie seen -- node ids + equiv members -- to resolve through
      the normalizer, whose cliques are the equivalence backbone; see stage 1b).

    Category strings are interned so the facts dict stays compact (there are only
    ~150 distinct Biolink categories).
    """
    inherited_cats: dict[str, set[str]] = {}
    node_taxon: dict[str, str] = {}
    node_ids: set[str] = set()
    # seeds: every id we've seen -> a bitmask of the sources that provided it (a node
    # id or an equiv-list member). Doubles as the normalizer seed set (its keys) and
    # the provenance for retaining EVERY id as a node (singleton if it never merges).
    seeds: dict[str, int] = {}
    # per-source infores provided_by (so a retained bare id can carry real provenance).
    source_provided_by: dict[str, set[str]] = defaultdict(set)

    def record_aliases(edge: dict, source: str, bit: int, ev) -> None:
        """Attach the ids an aggregator STARTED from that exist nowhere else in the graph.

        An aggregator edge stores "I resolved X to Y" (see ``original_alias_pairs``). When X appears in NO
        source's nodes or equivalency lists, that pairing is the only thing that can place it, so X is seeded
        like an equiv-list member (normalizer name/category/taxon, real provenance, a node of its own if it
        never merges) and the X/Y pairing is written as evidence weighted like a two-id list from that
        aggregator -- merge strength, so X joins Y's cluster. This is how robokop's gtex edges, which originate
        at ``HGVS:`` ids no node carries, reach their ``CAID:`` variants.

        When X IS already in the graph, the alias is skipped entirely: some source's list or a normalizer
        clique already places X, and an alias could only repeat that or override it. Every bad merge aliases
        caused came from that case -- thousands of ids fanning into one canonical, or originals paired with
        the wrong endpoint -- because only there can an alias join two entities that already exist.

        "Already in the graph" is decided after EVERY source's nodes have been read (the node pass runs to
        completion before any edges), so it doesn't depend on the order sources are processed in.
        """
        for original, canonical in original_alias_pairs(edge, source):
            mask = seeds.get(original, 0)
            if mask and not mask & ALIAS_SEED_BIT:
                continue  # already placed by some source's nodes or lists
            seeds[original] = mask | bit | ALIAS_SEED_BIT
            # The original inherits the canonical id's categories the same way an equiv-list member
            # does (the node pass). Without it, an id from a vocabulary the normalizer doesn't
            # know (HGVS, CAID) would fall through to NamedThing -- a guardrail wildcard -- so the
            # branch guardrail would go inert on exactly the ids this is introducing.
            canonical_cats = inherited_cats.get(canonical)
            if canonical_cats:
                inherited_cats.setdefault(original, set()).update(canonical_cats)
            if is_coarser_than_canonical(original, canonical):
                # e.g. a bare rsid (position) stored on one CAID (allele): seeded and typed above so the edge
                # can remap onto the position it was asserted about, but never merged into that one allele.
                continue
            ev_alias = alias_evidence(original, canonical, source, weights)
            if ev_alias is not None:
                ev.write(_evidence_row(*ev_alias, kind=f"alias:{source}"))

    sources = sorted(config.all_harmonized_paths_resolved.items())
    with open(evidence_path, "w") as ev, open(names_path, "w") as nm:
        # Pass 1: every source's NODES (ids, equivalency lists, categories, taxa, names).
        for source, (nodes_path, _edges_path) in sources:
            bit = 1 << source_bits[source]
            if Path(nodes_path).exists():
                for node in stream_nodes_from_jsonl(Path(nodes_path)):
                    node_id = node.get(NODE_ID)
                    if not node_id:
                        continue
                    node_ids.add(node_id)
                    seeds[node_id] = seeds.get(node_id, 0) | bit
                    source_provided_by[source].update(node.get(NODE_PROVIDED_BY) or ())
                    equiv_ids = node.get(NODE_EQUIVALENT_IDS) or []
                    for equiv_id in equiv_ids:
                        seeds[equiv_id] = seeds.get(equiv_id, 0) | bit
                    # Equivalency-clique evidence from EVERY source, weighted per source:
                    # native curated lists are strong (>=tau, merge on their own). The
                    # canonicalized aggregators are PREFIX-CAPPED -- ids of a prefix a list holds
                    # in bulk get only a weak link to the listing node, the rest merge (see
                    # ERWeights.max_ids_per_prefix) -- and share the "sri_nn_derived"
                    # source group so their Babel echo counts once. The clean, current
                    # cross-ontology backbone comes from the normalizer's cliques (stage 1b);
                    # aggregator lists recover the mappings it doesn't know.
                    for evidence in clique_evidence(equiv_ids, source, weights, head=node_id):
                        ev.write(_evidence_row(*evidence, kind=f"equiv:{source}"))
                    cats = node.get(NODE_CATEGORIES) or []
                    branches = families.branches(cats) if cats else ALL_FAMILIES
                    if cats and branches is not ALL_FAMILIES and len(branches) == 1:
                        interned = tuple(sys.intern(c) for c in cats)
                        inherited_cats.setdefault(node_id, set()).update(interned)
                        for equiv_id in equiv_ids:
                            inherited_cats.setdefault(equiv_id, set()).update(interned)
                    taxon = node.get(NODE_TAXON)
                    if taxon:
                        node_taxon[node_id] = sys.intern(taxon)
                    # Name rows for name-similarity are PER-ID. A NATIVE source names
                    # its own id (attributable). A CANONICALIZED aggregator node's name
                    # is the clique's PREFERRED label -- not reliably the canonical id's
                    # own name -- so we do NOT emit it; those ids are named per-id by
                    # their NN label instead (stage 1b).
                    if source not in CANONICALIZED_AGGREGATOR_SOURCES:
                        name_norm = normalize_name(node.get(NODE_NAME))
                        if not is_droppable(name_norm, min_length=weights.min_name_length, stoplist=DEFAULT_STOPLIST):
                            nm.write(f"{name_norm}{SEP}{node_id}\n")

        # Pass 2: every source's EDGES -- aliases (which need pass 1 complete) and match predicates.
        # Match-predicate (close/exact/same_as) edges are match-graph evidence, but
        # only on their ORIGINAL endpoints (see _original_endpoints): a Babel-
        # canonicalized endpoint would just re-import Babel's clustering. KG2 also
        # needs its close_match down-weighted where subclass edges co-occur, so it
        # gets a dedicated two-phase writer.
        for source, (_nodes_path, edges_path) in sources:
            bit = 1 << source_bits[source]
            if not Path(edges_path).exists():
                continue
            if source == "kg2":
                _write_kg2_match_evidence(Path(edges_path), weights, ev)
                continue
            takes_aliases = source in ALIAS_EVIDENCE_SOURCES
            for edge in stream_edges_from_jsonl(Path(edges_path)):
                if takes_aliases:
                    record_aliases(edge, source, bit, ev)
                predicate = edge.get(EDGE_PREDICATE, "")
                if weights.predicate_weight(predicate) is None:
                    continue  # not a usable match predicate; skip before un-canon work
                endpoints = _original_endpoints(edge, source)
                if endpoints is None:
                    continue  # canonicalized aggregator without a known un-canonicalizer
                primary_ks = _primary_ks(edge)
                for subject, object_ in endpoints:
                    ev_edge = match_predicate_evidence(
                        subject, object_, predicate, source, weights, primary_ks=primary_ks
                    )
                    if ev_edge is not None:
                        ev.write(_evidence_row(*ev_edge, kind=f"match:{source}"))
    return inherited_cats, node_taxon, node_ids, seeds, source_provided_by


def _stage1b_normalizer_evidence_and_names(
    nodenorm: NodeNormClient,
    seeds: set[str],
    evidence_path: Path,
    names_path: Path,
    weights: ERWeights,
) -> None:
    """Resolve every seed through the normalizer (harvesting whole cliques, so each
    clique costs one query, not one per member), then append two things:

    * **equivalence evidence** from the normalizer's cliques -- the clean, current
      cross-ontology backbone that replaces the aggregators' baked-in lists;
    * **per-id name rows** from each id's own NN label -- so name-similarity links
      INDIVIDUAL identifiers by their INDIVIDUAL names (not by some source node's
      primary name), which is the only correct granularity for it.
    """
    logging.info("entity_resolution: resolving %d seed curies via node normalizer (harvesting cliques)", len(seeds))
    nodenorm.resolve(seeds)
    n_cliques = 0
    with open(evidence_path, "a") as ev:
        for canonical, members in nodenorm.iter_cliques():
            for evidence in clique_evidence(members, "nn", weights, head=canonical):
                ev.write(_evidence_row(*evidence, kind="nn"))
            n_cliques += 1
    n_names = 0
    with open(names_path, "a") as nm:
        for curie, label in nodenorm.iter_labels():
            name_norm = normalize_name(label)
            if not is_droppable(name_norm, min_length=weights.min_name_length, stoplist=DEFAULT_STOPLIST):
                nm.write(f"{name_norm}{SEP}{curie}\n")
                n_names += 1
    logging.info(
        "entity_resolution: appended %d normalizer cliques + %d per-id name rows", n_cliques, n_names
    )


def _primary_ks(edge: dict) -> str | None:
    """An edge's primary knowledge source, as a single id (harmonized edges may carry a list)."""
    primary_ks = edge.get(EDGE_PRIMARY_KS)
    if isinstance(primary_ks, list):
        return primary_ks[0] if primary_ks else None
    return primary_ks or None


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

    KG2's originals contribute no alias evidence (see uncanonicalize.ALIAS_EVIDENCE_SOURCES).
    """
    hierarchical_counts: dict[tuple[str, str], int] = defaultdict(int)
    match_pairs: list[tuple[str, str, str, str | None]] = []  # (a, b, predicate, primary_ks), a <= b
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
                match_pairs.append((a, b, predicate, _primary_ks(edge)))

    for a, b, predicate, primary_ks in match_pairs:
        base = weights.predicate_weight(predicate)
        weight = _subclass_penalized_weight(base, hierarchical_counts.get((a, b), 0), weights.subclass_penalty_decay)
        # Per-primary-KS group, so parallel close_matches on this pair from different KSes sum.
        group = weights.predicate_group("kg2", primary_ks)
        ev.write(_evidence_row(a, b, group, weight, kind="match:kg2"))


def _stage1c_append_name_similarity(names_path: Path, evidence_path: Path, weights: ERWeights, temp_dir: Path) -> None:
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
                        ev.write(_evidence_row(a, b, NAME_SIMILARITY_GROUP, w, kind="name_sim"))

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
    """Combine evidence per CURIE pair (max within source group, sum across)
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
                a, b, group, weight_s = line.rstrip("\n").split(SEP)[:4]  # 5th column (evidence kind) is diagnostic
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
    inherited_cats: dict[str, set[str]],
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

    # One category SET per match-graph node (may hold several biolink categories —
    # we keep all of them, and the node's family is their union), from a chain that
    # never leaves it untyped:
    #   1. the node normalizer's per-id type LIST (source of truth; types each id
    #      individually, so a conflated MONDO node becomes just Disease, while a
    #      genuine bridge id keeps both types — e.g. [ChemicalEntity, Protein]);
    #   2. else the categories inherited from the single-family nodes that list this id
    #      in their equivalency list — the full inherited set, but ONLY if it resolves
    #      to a single family; cross-family disagreement (different aggregators typing
    #      it differently) is ambiguous and falls through;
    #   3. else the prefix->category backup (infer_category) — a LAST-resort guess that
    #      must never pre-empt a real source category, so it runs AFTER inheritance, not
    #      inside the normalizer client;
    #   4. else NamedThing (a guardrail wildcard).
    # This keeps the branch guardrail from going inert on the huge fraction of ids
    # that only ever appear as equiv-list members, without re-importing Babel's
    # mis-typing (multi-family nodes never propagate in stage 1).
    all_curies = [uniques[i] for i in range(num_nodes)]
    logging.info("entity_resolution: resolving %d match-graph node categories/taxa (cached)", len(all_curies))
    resolved = nodenorm.resolve(all_curies)  # cache hits after stage 1b; no new API calls
    mg_categories: dict[str, tuple[str, ...]] = {}
    mg_taxon: dict[str, str] = {}
    bare_names: dict[str, str] = {}
    for curie in all_curies:
        norm = resolved.get(curie)
        cats = tuple(norm.categories) if norm and norm.categories else ()  # 1. normalizer (source of truth)
        if not cats:
            inherited = inherited_cats.get(curie)  # 2. inherited from single-family source equiv-lists
            if inherited:
                branches = families.branches(tuple(inherited))
                if branches is not ALL_FAMILIES and len(branches) == 1:
                    cats = tuple(sorted(inherited))  # single-family inherited category
        if not cats:
            inferred = infer_category(curie)  # 3. prefix backup -- LAST resort, never before source inheritance
            if inferred:
                cats = (inferred,)
        if not cats:
            cats = ("biolink:NamedThing",)  # 4. untyped/ambiguous -> NamedThing (a guardrail wildcard)
        mg_categories[curie] = cats
        # Taxon precedence (mirrors category): 1. normalizer (source of truth);
        # 2. the harmonized node's (source) taxon; 3. the single-species prefix backup
        # -- LAST, never before source; else untaxoned (a guardrail wildcard).
        taxa = norm.taxa if norm else ()
        if taxa:
            mg_taxon[curie] = taxa[0]
        elif curie in node_taxon:
            mg_taxon[curie] = node_taxon[curie]
        else:
            inferred_taxon = infer_taxon(curie)
            if inferred_taxon:
                mg_taxon[curie] = inferred_taxon
        if norm and norm.label and curie not in node_ids:
            bare_names[curie] = norm.label

    def info_provider(curie: str) -> NodeInfo:
        # Stamp the node with its resolved branch-FAMILY set (category -> family done
        # once here), so the guardrails cluster on families directly.
        return NodeInfo(
            curie=curie,
            branches=families.branches(mg_categories.get(curie, ())),
            taxon=mg_taxon.get(curie),
        )

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
                e for e in comp_edges if not cluster_violations([e[0], e[1]], info, guardrail_config)
            ]
            # Label propagation on the pruned component (no resolution parameter —
            # LP merges what's connected; the guardrails do the splitting).
            raw_clusters = label_propagation(member_curies, comp_edges, seed=seed)
            # Guardrails as the backstop, split until valid (catches transitive
            # conflicts the pairwise prune can't see). LP has no resolution to raise,
            # so there's no clustering-based splitter — greedy_valid_partition repairs.
            adjacency = _adjacency(comp_edges)
            checked: list[list[str]] = []
            for cluster in raw_clusters:
                checked.extend(
                    enforce_cluster(cluster, info, guardrail_config, splitter=None, adjacency=adjacency)
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


# --------------------------------------------------------------------------------------
# Stage 4: materialize canonical nodes (external group-by) -> node_id -> representative
# --------------------------------------------------------------------------------------


def _stage4_materialize(
    config,
    curie_to_cluster: dict[str, int],
    mg_categories: dict[str, tuple[str, ...]],
    node_ids: set[str],
    seeds: dict[str, int],
    source_provided_by: dict[str, set[str]],
    source_bits: dict[str, int],
    inherited_cats: dict[str, set[str]],
    nodenorm: NodeNormClient,
    ranking: PrefixRanking,
    families: BranchFamilies,
    biolink,
    temp_dir: Path,
) -> dict[str, str]:
    """Stream harmonized nodes, group by cluster on disk, reconcile one cluster at a
    time, and write the canonical nodes file. Returns ``node_id -> representative``.

    EVERY id a source provided is materialized -- merged into its cluster, or emitted
    as its own SINGLETON if it never merged. Bare ids (equiv-list members with no
    harmonized node) become synthetic member dicts carrying their normalizer
    name/category/taxon and real provenance (decoded from the ``seeds`` bitmask), so
    nothing a source provided is ever dropped. Peak memory is one cluster's members.
    """
    keyed = temp_dir / "er_s4_nodes_keyed.tmp"
    keyed_sorted = temp_dir / "er_s4_nodes_keyed_sorted.tmp"
    node_id_to_rep: dict[str, str] = {}
    bit_to_source = {bit: src for src, bit in source_bits.items()}

    def provenance(mask: int) -> list[str]:
        provided: set[str] = set()
        for bit, src in bit_to_source.items():
            if mask & (1 << bit):
                provided.update(source_provided_by.get(src, ()) or {src})
        return sorted(provided)

    def label_of(curie: str) -> str | None:
        info = nodenorm.get(curie)
        return info.label if info else None

    def taxon_of(curie: str) -> str | None:
        info = nodenorm.get(curie)
        return info.taxa[0] if info and info.taxa else None

    def bare_category(curie: str) -> list[str]:
        cats = mg_categories.get(curie)  # computed in stage 3 for match-graph ids
        if cats:
            return sorted(cats)
        info = nodenorm.get(curie)  # isolated id: NN -> inherited -> prefix backup -> NamedThing
        if info and info.categories:
            return sorted(info.categories)
        inherited = inherited_cats.get(curie)
        if inherited:
            branches = families.branches(tuple(inherited))
            if branches is not ALL_FAMILIES and len(branches) == 1:
                return sorted(inherited)
        inferred = infer_category(curie)  # prefix backup only after source inheritance (see mg_categories)
        if inferred:
            return [inferred]
        return ["biolink:NamedThing"]

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
                    # Taxon: NN wins; else keep the source taxon already on the dict; else
                    # the single-species prefix backup (last resort, never before source).
                    nn_taxon = taxon_of(node_id)
                    if nn_taxon:
                        node[NODE_TAXON] = nn_taxon
                    elif not node.get(NODE_TAXON):
                        inferred_taxon = infer_taxon(node_id)
                        if inferred_taxon:
                            node[NODE_TAXON] = inferred_taxon
                    cid = curie_to_cluster.get(node_id)
                    key = f"c{cid}" if cid is not None else f"s:{node_id}"
                    writer.write([key, node])
            # Synthetic member dicts for EVERY bare id (equiv-list member with no
            # harmonized node): merged into its cluster if it merged, else its own
            # SINGLETON -- so every id a source provided survives to the graph.
            for curie, mask in seeds.items():
                if curie in node_ids:
                    continue  # has a harmonized node (already written above)
                cid = curie_to_cluster.get(curie)
                key = f"c{cid}" if cid is not None else f"s:{curie}"
                synthetic = {
                    NODE_ID: curie,
                    NODE_CATEGORIES: bare_category(curie),
                    NODE_PROVIDED_BY: provenance(mask),
                }
                info = nodenorm.get(curie)
                if info and info.label:
                    synthetic[NODE_NAME] = info.label
                if info and info.taxa:
                    synthetic[NODE_TAXON] = info.taxa[0]
                else:
                    inferred_taxon = infer_taxon(curie)  # bare id: no source taxon; prefix backup last
                    if inferred_taxon:
                        synthetic[NODE_TAXON] = inferred_taxon
                writer.write([key, synthetic])
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
                node = materialize_cluster(group, ranking, families, label_of=label_of)
                # No drop rule: every id a source provided is kept (bare ids now carry
                # real provenance, so nothing is provenance-less). A bare-only cluster
                # becomes its own node instead of vanishing.
                if biolink is not None and node.get(NODE_CATEGORIES):
                    leaves = biolink.filter_to_leaf_categories(node[NODE_CATEGORIES])
                    if leaves:
                        node[NODE_CATEGORIES] = sorted(leaves)
                # Every node needs a category; an untyped survivor becomes NamedThing.
                if not node.get(NODE_CATEGORIES):
                    node[NODE_CATEGORIES] = ["biolink:NamedThing"]
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


def _stage_banner(msg: str) -> None:
    """Log a prominent, greppable banner so the numbered ER stages stand out in the very
    long build log (otherwise they're lost among hundreds of thousands of INFO lines)."""
    logging.info("#" * 100)
    logging.info("###  %s", msg)


# Every cluster at least this large gets a breakdown of the evidence that built it, written to
# ``<integrated debug dir>/oversized_clusters.jsonl``. 2.1.1's largest was 1,455 members, so a healthy build writes a
# handful; a regression (the 56,515-member clusters of the first size-aware build) shows what caused it.
OVERSIZED_EVIDENCE_REPORT_MIN_SIZE = 1000
OVERSIZED_REPORT_FILENAME = "oversized_clusters.jsonl"
OVERSIZED_REPORT_SAMPLE_IDS = 25


def _report_oversized_cluster_evidence(
    evidence_path: Path, curie_to_cluster: dict[str, int], report_path: Path
) -> None:
    """Write, for every cluster of OVERSIZED_EVIDENCE_REPORT_MIN_SIZE+ members, how many evidence rows of each kind link
    members INSIDE it, plus its id-prefix mix and a sample of ids -- one JSON line per cluster, largest first.

    One pass over the evidence file, and only when such a cluster exists. The kind column (``equiv:<source>``,
    ``alias:<source>``, ``match:<source>``, ``nn``, ``name_sim``) is what separates a normalizer clique from an
    aggregator list; the de-correlation group alone can't.
    """
    sizes = Counter(curie_to_cluster.values())
    ranked = [cid for cid, n in sizes.most_common() if n >= OVERSIZED_EVIDENCE_REPORT_MIN_SIZE]
    if not ranked:
        return
    wanted = set(ranked)
    kinds: dict[int, Counter] = {cid: Counter() for cid in ranked}
    prefixes: dict[int, Counter] = {cid: Counter() for cid in ranked}
    samples: dict[int, list[str]] = defaultdict(list)
    for curie, cid in curie_to_cluster.items():
        if cid in wanted:
            prefixes[cid][curie.split(":", 1)[0]] += 1
            if len(samples[cid]) < OVERSIZED_REPORT_SAMPLE_IDS:
                samples[cid].append(curie)
    with open(evidence_path) as evidence:
        for line in evidence:
            a, b, _group, _weight, kind = line.rstrip("\n").split(SEP)
            cid = curie_to_cluster.get(a)
            if cid in wanted and curie_to_cluster.get(b) == cid:
                kinds[cid][kind] += 1
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with jsonlines.open(report_path, "w") as report:
        for cid in ranked:
            report.write(
                {
                    "members": sizes[cid],
                    "evidence_rows_inside_by_kind": dict(kinds[cid].most_common()),
                    "id_prefixes": dict(prefixes[cid].most_common()),
                    "sample_ids": sorted(samples[cid]),
                }
            )
    logging.warning(
        "entity_resolution: %d cluster(s) have %d+ members (largest %d); evidence breakdown written to %s",
        len(ranked),
        OVERSIZED_EVIDENCE_REPORT_MIN_SIZE,
        sizes[ranked[0]],
        report_path,
    )


def _debug_dir(config) -> Path:
    """Where diagnostics go: the build's integrated debug dir (falls back to the integrated dir for ad-hoc configs)."""
    return Path(getattr(config, "integrated_debug_dir", None) or config.integrated_dir)


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

    # Deterministic source -> bit index, for encoding per-id provenance in ``seeds``.
    source_bits = {src: i for i, src in enumerate(sorted(config.all_harmonized_paths_resolved))}

    config.er_nodenorm_cache_path.parent.mkdir(parents=True, exist_ok=True)
    nodenorm = NodeNormClient(config.er_nodenorm_cache_path)
    try:
        t = time.perf_counter()
        _stage_banner("ER STAGE 1 -- streaming harmonized nodes/edges -> match evidence, names, guardrail facts")
        inherited_cats, node_taxon, node_ids, seeds, source_provided_by = _stage1_write_evidence_and_facts(
            config, weights, families, evidence_path, names_path, source_bits
        )
        _stage_banner(
            f"ER STAGE 1 DONE ({time.perf_counter() - t:.1f}s) -- {len(node_ids)} harmonized node ids, "
            f"{len(seeds)} total ids (incl. equiv-list members)"
        )

        # 1b (normalizer cliques + per-id NN name rows) runs before 1c so name-similarity
        # groups over EVERY id's individual name, not just harmonized primary names.
        t = time.perf_counter()
        _stage_banner("ER STAGE 1b -- fetching Node Normalizer cliques + per-id names (the equivalence backbone)")
        _stage1b_normalizer_evidence_and_names(nodenorm, seeds, evidence_path, names_path, weights)
        _stage_banner(f"ER STAGE 1b DONE ({time.perf_counter() - t:.1f}s)")

        t = time.perf_counter()
        _stage_banner("ER STAGE 1c -- adding name-similarity match evidence")
        _stage1c_append_name_similarity(names_path, evidence_path, weights, temp_dir)
        _stage_banner(f"ER STAGE 1c DONE ({time.perf_counter() - t:.1f}s)")

        t = time.perf_counter()
        _stage_banner("ER STAGE 2 -- accumulating + tau-filtering weighted match pairs")
        n_pairs = _stage2_accumulate_pairs(evidence_path, pairs_path, weights, temp_dir)
        _stage_banner(f"ER STAGE 2 DONE ({time.perf_counter() - t:.1f}s) -- {n_pairs} pairs above tau")

        t = time.perf_counter()
        _stage_banner("ER STAGE 3 -- clustering (connected components -> label propagation -> guardrails)")
        curie_to_cluster, histogram, _bare_names, mg_categories = _stage3_cluster(
            pairs_path,
            weights,
            families,
            guardrail_config,
            inherited_cats,
            node_taxon,
            node_ids,
            nodenorm,
            DEFAULT_SEED,
        )
        _stage_banner(
            f"ER STAGE 3 DONE ({time.perf_counter() - t:.1f}s) -- "
            f"{len(curie_to_cluster)} ids -> {len(set(curie_to_cluster.values()))} clusters"
        )
        _log_histogram(histogram)
        _report_oversized_cluster_evidence(
            evidence_path, curie_to_cluster, _debug_dir(config) / OVERSIZED_REPORT_FILENAME
        )

        t = time.perf_counter()
        _stage_banner(f"ER STAGE 4 -- materializing canonical nodes -> {config.integrated_nodes_path}")
        node_id_to_rep = _stage4_materialize(
            config,
            curie_to_cluster,
            mg_categories,
            node_ids,
            seeds,
            source_provided_by,
            source_bits,
            inherited_cats,
            nodenorm,
            ranking,
            families,
            biolink,
            temp_dir,
        )
        _stage_banner(
            f"ER STAGE 4 DONE ({time.perf_counter() - t:.1f}s) -- "
            f"{len(node_id_to_rep)} node ids -> {len(set(node_id_to_rep.values()))} canonical nodes"
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
