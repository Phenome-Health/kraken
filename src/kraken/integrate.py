"""Graph integration: clustering-based entity resolution + edge merging.

Node resolution is handled by ``kraken.entity_resolution`` (see that package and
``docs/entity_resolution_plan.md``): it clusters the match graph and writes the
canonical nodes file, returning a ``node_id -> representative_curie`` map. Edges
are then resolved through that map and merged by the existing order-independent
external sort. The legacy ``primary_source`` / ``can_merge_existing_nodes`` node
merge has been fully replaced.
"""

import json
import logging
import os
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import jsonlines

from kraken.biolink_client import BiolinkClient
from kraken.config import KrakenConfig
from kraken.entity_resolution.build import resolve_entities
from kraken.entity_resolution.uncanonicalize import (
    kg2_pre_id_triples,
    original_endpoints,
)
from kraken.schema import EdgeModel
from kraken.utils.constants import (
    CROSS_CLUSTER_EQUIVALENCE_PREDICATE,
    EDGE_AGENT_TYPE,
    EDGE_AGGREGATOR_KS,
    EDGE_ATTRIBUTES,
    EDGE_KNOWLEDGE_LEVEL,
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_PRIMARY_KS,
    EDGE_SUBJECT,
    KNOWLEDGE_ASSERTION,
    KRAKEN_SOURCE_ID,
    NODE_EQUIVALENT_IDS,
    NODE_ID,
    NODE_PROVIDED_BY,
    NOT_PROVIDED,
    SAME_AS_PREDICATE,
)
from kraken.utils.general import create_edge_key, to_list
from kraken.utils.kg_io import remove_file, stream_edges_from_jsonl, stream_nodes_from_jsonl


def integrate_sources(config: KrakenConfig, biolink: BiolinkClient):
    """Resolve entities into canonical nodes, then merge edges across all sources."""
    # Babel is the equivalence backbone and the per-id name/category/taxon authority (see entity_resolution.build);
    # without it entity resolution would still run, but quietly lose most cross-ontology merges.
    if "babel" not in config.sources_to_use:
        raise ValueError(
            "Integration requires the 'babel' source: entity resolution takes its cliques and per-id names, "
            "categories and taxa from it. Include babel (and harmonize it) before integrating."
        )

    config.integrated_dir.mkdir(parents=True, exist_ok=True)
    config.integrated_debug_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Starting source integration (clustering-based entity resolution)...")

    # Phase 1: cluster the match graph -> canonical nodes file + node_id -> representative map.
    node_id_to_rep = resolve_entities(config, biolink)
    assert node_id_to_rep, "entity resolution produced no nodes"

    # Phase 2: resolve edge endpoints through that map and merge duplicate edges.
    integrate_edges(node_id_to_rep, config)

    logging.info(f"Integration complete! Unified KG saved to {config.integrated_dir}")


# Field separator for the temporary key-sorted edge file. Safe because json.dumps escapes any tabs or
# newlines inside string values, so neither byte ever appears within a serialized edge.
_EDGE_SORT_SEP = "\t"


def integrate_edges(node_map: dict[str, str], config: KrakenConfig):
    """Merge edges across ALL sources using a disk-based external sort.

    Edges sharing an edge key are merged into a single edge. Because aggregator_knowledge_source is not
    part of the key, the same assertion arriving via different aggregators (e.g. kg2 and robokop) collapses
    into one edge, with its aggregator_knowledge_source list union-merged.

    Rather than holding every mergeable edge in memory, we stream all edges to a temp file keyed by edge
    key, sort that file on disk, then merge each run of same-key edges in a single streaming pass. Peak
    memory is therefore just one group of same-key edges, regardless of graph size.
    """
    assert node_map

    keyed_edges_path = config.integrated_dir / "edges_keyed.tmp.tsv"
    sorted_edges_path = config.integrated_dir / "edges_keyed_sorted.tmp.tsv"
    try:
        _write_keyed_edges(node_map, config, keyed_edges_path)
        _sort_file_by_key(keyed_edges_path, sorted_edges_path, temp_dir=config.integrated_dir)
        total_edges, num_merged = _merge_sorted_edges(sorted_edges_path, config)
        logging.info(f"Wrote {total_edges} integrated edges ({num_merged} merged from multiple source edges)")
    finally:
        remove_file(keyed_edges_path)
        remove_file(sorted_edges_path)


def _write_keyed_edges(node_map: dict[str, str], config: KrakenConfig, keyed_edges_path: Path):
    """Stream every source's edges to a temp file as '<edge_key>\\t<edge_json>' lines, resolving each
    edge's subject/object to canonical (representative) IDs first so the keys reflect the integrated graph.

    Edges are remapped by their ORIGINAL endpoints (un-canonicalizing aggregators): a canonicalizing
    aggregator stores Babel-canonical endpoints, but our clustering diverges from Babel, so mapping a
    canonical endpoint would attach the edge to the wrong node. One stored edge can carry several
    original subject/object pairs (e.g. kg2), so it fans out into several edges. Sources we can't
    un-canonicalize fall back to their stored endpoints.

    KG2's originals are NOT assumed to be listed in the stored edge's direction -- KG2 re-orients edges when it
    normalizes an inverse relation, so ~12% of its original pairs run backwards (see ``_orient``). Other sources'
    originals are taken in the order they're recorded."""
    orphaned = self_loops = 0
    orientation = Counter()
    votes: dict[tuple[str, str], Counter] = defaultdict(Counter)  # (predicate, original relation) -> orientations
    spool_path = keyed_edges_path.with_suffix(".undetermined.tmp")

    def write(edge: dict, rep_subj: str, rep_obj: str) -> None:
        nonlocal self_loops
        if rep_subj == rep_obj:
            self_loops += 1  # endpoints merged into one node -> self-loop, drop
            return
        resolved = {**edge, EDGE_SUBJECT: rep_subj, EDGE_OBJECT: rep_obj}
        keyed_file.write(f"{create_edge_key(resolved)}{_EDGE_SORT_SEP}{json.dumps(resolved)}\n")

    try:
        with open(keyed_edges_path, "w") as keyed_file, open(spool_path, "w") as spool:
            for source_name in config.sources_to_use:
                logging.info(f"Writing keyed edges from {source_name}..")
                _, edges_file = config.all_harmonized_paths_resolved[source_name]
                for edge in stream_edges_from_jsonl(edges_file):
                    if source_name == "babel" and edge.get(EDGE_PREDICATE) == SAME_AS_PREDICATE:
                        # A Babel clique edge. Within one cluster it becomes a self-loop and is dropped below; one
                        # that survives joins ids entity resolution kept apart, so it can't claim same_as.
                        edge = {**edge, EDGE_PREDICATE: CROSS_CLUSTER_EQUIVALENCE_PREDICATE}
                    # KG2's originals carry their relation (needed to settle undetermined orientations)
                    triples = kg2_pre_id_triples(edge) if source_name == "kg2" else []
                    if not triples:
                        pairs = original_endpoints(edge, source_name)
                        if pairs is None:  # canonicalized aggregator with no recoverable originals -> stored endpoints
                            pairs = [(edge.get(EDGE_SUBJECT), edge.get(EDGE_OBJECT))]
                        triples = [(subj, None, obj) for subj, obj in pairs]
                    for subj_id, relation, obj_id in triples:
                        rep_a, rep_b = node_map.get(subj_id), node_map.get(obj_id)
                        if rep_a is None or rep_b is None:
                            orphaned += 1  # an endpoint that never became a node -> skip this edge
                            continue
                        if source_name != "kg2":
                            # Only KG2 re-orients edges relative to its recorded originals. ROBOKOP and
                            # Translator write original_subject/original_object to match their own fields,
                            # and a native source's endpoints are its own -- orienting those could only flip
                            # a correct edge when an original happens to be clustered with the other end.
                            write(edge, rep_a, rep_b)
                            continue
                        direction = _orient(rep_a, rep_b, edge, node_map)
                        orientation[direction] += 1
                        if direction == "undetermined":
                            if relation is not None:  # decide after the pass, from how this relation oriented
                                spool.write(json.dumps([edge, rep_a, rep_b, relation]) + "\n")
                            else:
                                write(edge, rep_a, rep_b)
                            continue
                        if relation is not None:
                            votes[(edge.get(EDGE_PREDICATE), relation)][direction] += 1
                        if direction == "swapped":
                            write(edge, rep_b, rep_a)
                        else:
                            write(edge, rep_a, rep_b)

            spool.close()
            by_relation = Counter()
            with open(spool_path) as undetermined:
                for line in undetermined:
                    edge, rep_a, rep_b, relation = json.loads(line)
                    tally = votes.get((edge.get(EDGE_PREDICATE), relation))
                    if tally and tally["swapped"] > tally["aligned"]:
                        by_relation["swapped"] += 1
                        write(edge, rep_b, rep_a)
                    else:
                        by_relation["aligned" if tally else "no evidence (kept original order)"] += 1
                        write(edge, rep_a, rep_b)
            _write_equivalence_edges(node_map, config, keyed_file)
    finally:
        remove_file(spool_path)

    if orientation:
        logging.info(
            "Oriented KG2 original endpoints by which stored endpoint's cluster each lands in: "
            "%d aligned, %d swapped (re-oriented), %d undetermined -> resolved by relation: %s",
            orientation["aligned"],
            orientation["swapped"],
            orientation["undetermined"],
            dict(by_relation),
        )
    if orphaned:
        logging.warning("Skipped %d edge endpoints with no node mapping (orphans)", orphaned)
    if self_loops:
        logging.info("Dropped %d self-loop edges (endpoints merged into one node)", self_loops)


def _orient(rep_a: str, rep_b: str, edge: dict, node_map: dict[str, str]) -> str:
    """Which way an original pair runs relative to its stored edge: "aligned", "swapped", or "undetermined".

    Decided by where the originals landed: an original in the stored SUBJECT's cluster is the subject. Clusters
    rather than id equality, because the originals are different ids from the (Babel-canonical) stored endpoints.
    Undetermined when neither original shares a cluster with either stored endpoint, or when the evidence points
    both ways (e.g. the stored endpoints themselves share a cluster).
    """
    stored_subj, stored_obj = node_map.get(edge.get(EDGE_SUBJECT)), node_map.get(edge.get(EDGE_OBJECT))
    aligned = rep_a == stored_subj or rep_b == stored_obj
    swapped = rep_a == stored_obj or rep_b == stored_subj
    if aligned and not swapped:
        return "aligned"
    if swapped and not aligned:
        return "swapped"
    return "undetermined"


def _cross_cluster_edge(subject: str, object_: str, primary_ks: str, aggregator_ks: list[str]) -> dict:
    """A close_match edge between two clusters asserted equivalent but kept apart (symmetric: subject/object
    ordered so A~B and B~A collapse). See CROSS_CLUSTER_EQUIVALENCE_PREDICATE for why not same_as.

    knowledge_level = knowledge_assertion (an asserted equivalence, not a prediction/
    statistic). agent_type = not_provided: these edges are synthesized from equiv-list
    co-membership, so the agent that originally asserted the equivalence is unknown (it
    varies by source), and we don't claim one. KRAKEN is recorded as the last aggregator in the chain, as it is
    on every directly ingested edge (see BaseHarmonizer.create_edge) -- the source asserted the equivalence, but
    the EDGE only exists because our entity resolution kept the two ids apart."""
    subject, object_ = sorted((subject, object_))
    edge = {
        EDGE_SUBJECT: subject,
        EDGE_OBJECT: object_,
        EDGE_PREDICATE: CROSS_CLUSTER_EQUIVALENCE_PREDICATE,
        EDGE_PRIMARY_KS: primary_ks,
        EDGE_KNOWLEDGE_LEVEL: KNOWLEDGE_ASSERTION,
        EDGE_AGENT_TYPE: NOT_PROVIDED,
    }
    edge[EDGE_AGGREGATOR_KS] = list(dict.fromkeys([*aggregator_ks, KRAKEN_SOURCE_ID]))
    return edge


def _write_equivalence_edges(node_map: dict[str, str], config: KrakenConfig, keyed_file):
    """Retain the equivalence signal as real edges: for every asserted equivalence whose
    two ids ended up in DIFFERENT clusters, emit a ``close_match`` edge between their
    representatives (same-cluster assertions collapse to self-loops and are dropped).
    So e.g. TP53 protein-isoforms that don't merge into the main TP53 node stay LINKED
    to it. These come from each source's equiv-lists (primary KS = that source); Babel's cliques are ordinary
    same_as edges, turned into close_match in _write_keyed_edges."""
    written = 0

    def emit(rep_a: str, rep_b: str, primary_ks: str, aggregator_ks: list[str]) -> int:
        if rep_a == rep_b:  # same cluster -> self-loop, skip
            return 0
        edge = _cross_cluster_edge(rep_a, rep_b, primary_ks, aggregator_ks)
        keyed_file.write(f"{create_edge_key(edge)}{_EDGE_SORT_SEP}{json.dumps(edge)}\n")
        return 1

    # (1) each source's equiv lists (the source asserts node_id ~ each of its members)
    for source_name in config.sources_to_use:
        if source_name == "babel":
            continue  # one node per id, so no lists: its equivalences are edges (see _write_keyed_edges)
        nodes_file, _ = config.all_harmonized_paths_resolved[source_name]
        for node in stream_nodes_from_jsonl(nodes_file):
            node_id = node.get(NODE_ID)
            rep_a = node_map.get(node_id)
            if rep_a is None:
                continue
            provided = node.get(NODE_PROVIDED_BY) or [source_name]
            primary_ks, aggregator_ks = provided[0], list(provided[1:])
            for member in node.get(NODE_EQUIVALENT_IDS) or []:
                rep_b = node_map.get(member)
                if rep_b is not None:
                    written += emit(rep_a, rep_b, primary_ks, aggregator_ks)

    logging.info("Wrote %d cross-cluster close_match equivalence edges", written)


def _sort_file_by_key(input_path: Path, output_path: Path, temp_dir: Path):
    """Externally sort a '<key>\\t<json>' file by its key column, using bounded memory. Byte ordering
    (LC_ALL=C) keeps it deterministic; -T keeps the sort's spill files on our (large) output volume."""
    logging.info("Sorting keyed edges on disk (external sort)..")
    subprocess.run(
        ["sort", "-t", _EDGE_SORT_SEP, "-k1,1", "-T", str(temp_dir), "-o", str(output_path), str(input_path)],
        check=True,
        env={**os.environ, "LC_ALL": "C"},
    )


def _merge_sorted_edges(sorted_edges_path: Path, config: KrakenConfig) -> tuple[int, int]:
    """Stream the key-sorted edges, merging each run of consecutive same-key edges into one. Merged edges
    are also written to a debug log. Returns (total_edges_written, num_merged_groups). Peak memory is a
    single same-key group."""
    total_written = 0
    num_merged = 0
    mergers_log = config.integrated_debug_dir / "edge_mergers.jsonl"
    with (
        open(sorted_edges_path) as sorted_file,
        jsonlines.open(config.integrated_edges_path, "w") as writer,
        jsonlines.open(mergers_log, "w") as mergers_writer,
    ):
        current_key = None
        group: list[dict] = []
        for line in sorted_file:
            key, _, edge_json = line.rstrip("\n").partition(_EDGE_SORT_SEP)
            if key != current_key and group:
                num_merged += _write_merged_group(group, writer, mergers_writer)
                total_written += 1
                group = []
            current_key = key
            group.append(json.loads(edge_json))
        if group:  # flush the final group
            num_merged += _write_merged_group(group, writer, mergers_writer)
            total_written += 1
    return total_written, num_merged


def _write_merged_group(group: list[dict], writer, mergers_writer) -> int:
    """Merge a group of same-key edges into a single edge and write it. Returns 1 if the group actually
    required merging (had more than one edge), else 0."""
    merged_edge = merge_edges(group)
    writer.write(merged_edge)
    if len(group) > 1:
        mergers_writer.write(merged_edge)
        return 1
    return 0


# Merging a group of same-key edges is done in ONE pass over the group, accumulating each property as it goes.
# It used to fold edges in one at a time, rebuilding every list-valued property from scratch on each fold -- O(n^2)
# in the group's size. That is invisible at the handful of edges a key normally has, but a single 417,750-edge group
# (thousands of wrongly-merged pathways, all pointing at the same object) needed ~87 billion set insertions and
# stalled integration for hours. The result is unchanged, except that merged lists now keep first-seen order
# instead of the arbitrary order a set gave them.


class _Union:
    """The distinct values of a merged property, accumulated in first-seen order in amortized O(1) per value.

    Unhashable values can't be de-duplicated, so they are all kept (as the old list merge did)."""

    __slots__ = ("_seen", "items")

    def __init__(self) -> None:
        self._seen: set = set()
        self.items: list = []

    def add(self, value: Any) -> None:
        for item in to_list(value):
            try:
                if item in self._seen:
                    continue
                self._seen.add(item)
            except TypeError:
                pass
            self.items.append(item)


class _DictFold:
    """A dict-valued property merged key by key (the one level of recursion the merge allows)."""

    __slots__ = ("values",)

    def __init__(self, first: dict) -> None:
        self.values: dict[str, Any] = dict(first)

    def add(self, value: dict) -> None:
        for key, item in value.items():
            self.values[key] = _fold(self.values.get(key), item, recursion_allowed=False, combine_flat_types=True)


def _fold(accumulated: Any, value: Any, *, recursion_allowed: bool = True, combine_flat_types: bool = False) -> Any:
    """Fold one more edge's ``value`` into a property's running ``accumulated`` value. Returns the new accumulation.

    The rules are the merge's long-standing ones: a missing value never replaces a present one; two dicts are merged
    key by key (top level only, where every entry combines); anything list-like -- or any flat value inside such a
    dict -- becomes the union of all values seen; otherwise the first value wins.
    """
    if value is None:
        return accumulated
    if accumulated is None:
        return value
    if isinstance(accumulated, _Union):
        accumulated.add(value)
        return accumulated
    if isinstance(accumulated, _DictFold):
        if isinstance(value, dict):
            accumulated.add(value)
            return accumulated
        accumulated = _materialize(accumulated)
    if isinstance(accumulated, dict) and isinstance(value, dict) and recursion_allowed:
        fold = _DictFold(accumulated)
        fold.add(value)
        return fold
    if (
        isinstance(accumulated, (set, list, dict, tuple))
        or isinstance(value, (set, list, dict, tuple))
        or combine_flat_types
    ):
        union = _Union()
        union.add(accumulated)
        union.add(value)
        return union
    return accumulated  # first flat value wins


def _materialize(accumulated: Any) -> Any:
    if isinstance(accumulated, _Union):
        return accumulated.items
    if isinstance(accumulated, _DictFold):
        return {key: _materialize(value) for key, value in accumulated.values.items()}
    return accumulated


_EDGE_KEY_PROPERTIES = EdgeModel.key_properties()


def merge_edges(group: list[dict]) -> dict:
    """Merge edges that share an edge key into one, in a single pass (linear in the group's total size).

    Key properties are identical across the group by definition, so the first edge's are kept. knowledge_level and
    agent_type take the first value that isn't not_provided. Attributes merge per source slot, and each slot's
    entries union. Every other property follows ``_fold``.
    """
    first = group[0]
    if len(group) == 1:
        return first
    accumulated: dict[str, Any] = {}
    for edge in group:
        for name, value in edge.items():
            if name in _EDGE_KEY_PROPERTIES:
                accumulated.setdefault(name, value)
            elif name in (EDGE_KNOWLEDGE_LEVEL, EDGE_AGENT_TYPE):
                if accumulated.get(name, NOT_PROVIDED) == NOT_PROVIDED:
                    accumulated[name] = value
            elif name == EDGE_ATTRIBUTES:
                slots = accumulated.setdefault(name, {})
                for source_slot, slot_value in (value or {}).items():
                    slots[source_slot] = _fold(slots.get(source_slot), slot_value)
            else:
                accumulated[name] = _fold(accumulated.get(name), value)
    merged = {}
    for name, value in accumulated.items():
        if name == EDGE_ATTRIBUTES:  # a plain dict of per-source accumulators, so materialize one level down
            merged[name] = {source_slot: _materialize(slot_value) for source_slot, slot_value in value.items()}
        else:
            merged[name] = _materialize(value)
    return merged
