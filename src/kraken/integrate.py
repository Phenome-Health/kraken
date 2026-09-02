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
from pathlib import Path
from typing import Any

import jsonlines

from kraken.biolink_client import BiolinkClient
from kraken.config import KrakenConfig
from kraken.entity_resolution.build import resolve_entities
from kraken.schema import EdgeModel
from kraken.utils.constants import (
    EDGE_AGENT_TYPE,
    EDGE_ATTRIBUTES,
    EDGE_KNOWLEDGE_LEVEL,
    EDGE_OBJECT,
    EDGE_SUBJECT,
    NODE_ATTRIBUTES,
    NODE_ID,
    NOT_PROVIDED,
)
from kraken.utils.general import create_edge_key, to_list
from kraken.utils.kg_io import remove_file, stream_edges_from_jsonl


def integrate_sources(config: KrakenConfig, biolink: BiolinkClient):
    """Resolve entities into canonical nodes, then merge edges across all sources."""
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
    edge's subject/object to canonical (representative) IDs first so the keys reflect the integrated graph."""
    with open(keyed_edges_path, "w") as keyed_file:
        for source_name in config.sources_to_use:
            logging.info(f"Writing keyed edges from {source_name}..")
            _, edges_file = config.all_harmonized_paths_resolved[source_name]
            for edge in stream_edges_from_jsonl(edges_file):
                resolve_to_canonical(edge, node_map)
                keyed_file.write(f"{create_edge_key(edge)}{_EDGE_SORT_SEP}{json.dumps(edge)}\n")


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
    merged_edge = group[0]
    for other_edge in group[1:]:
        merge_into_existing_edge(other_edge, merged_edge)
    writer.write(merged_edge)
    if len(group) > 1:
        mergers_writer.write(merged_edge)
        return 1
    return 0


def merge_into_existing_edge(new_edge: dict, existing_edge: dict):
    # NOTE: If edges are being merged, they must match on all properties included in the edge key

    # Merge knowledge_level, favoring values that aren't not_provided
    if existing_edge[EDGE_KNOWLEDGE_LEVEL] == NOT_PROVIDED:
        existing_edge[EDGE_KNOWLEDGE_LEVEL] = new_edge[EDGE_KNOWLEDGE_LEVEL]

    # Merge agent_type, favoring values that aren't not_provided
    if existing_edge[EDGE_AGENT_TYPE] == NOT_PROVIDED:
        existing_edge[EDGE_AGENT_TYPE] = new_edge[EDGE_AGENT_TYPE]

    # Merge any other properties as applicable (note: props included in edge key must be identical)
    for property_name, value in new_edge.items():
        if property_name not in EdgeModel.key_properties() | {EDGE_KNOWLEDGE_LEVEL, EDGE_AGENT_TYPE}:
            merge_property_into_existing(new_edge, existing_edge, property_name)


def resolve_to_canonical(edge: dict, node_map: dict[str, str]):
    subj_id = edge[EDGE_SUBJECT]
    obj_id = edge[EDGE_OBJECT]
    if subj_id in node_map and obj_id in node_map:
        edge[EDGE_SUBJECT] = node_map[subj_id]
        edge[EDGE_OBJECT] = node_map[obj_id]
    else:
        logging.warning(f"Skipping orphan edge: Edge between {subj_id} and {obj_id} is missing node mappings")


def merge_two_lists(list_a: list, list_b: list) -> list[Any]:
    # Merges two lists, retaining distinct values if hashable or otherwise just concatenating
    try:
        return list(set(list_a) | set(list_b))
    except Exception:
        return list_a + list_b


def merge_two_values(
    value_a: Any, value_b: Any, recursion_allowed: bool = True, combine_flat_types: bool = False
) -> Any:
    if value_a is None:
        return value_b
    elif value_b is None:
        return value_a
    elif isinstance(value_a, dict) and isinstance(value_b, dict) and recursion_allowed:
        # We recurse only on the top-level entries (no recursing beyond that, even if value is a dict)
        prop_names = set(value_a.keys()) | set(value_b.keys())
        merged_value = {
            prop_name: merge_two_values(
                value_a.get(prop_name), value_b.get(prop_name), recursion_allowed=False, combine_flat_types=True
            )
            for prop_name in prop_names
        }
        return merged_value
    elif (
        isinstance(value_a, (set, list, dict, tuple))
        or isinstance(value_b, (set, list, dict, tuple))
        or combine_flat_types
    ):
        value_a_list = to_list(value_a)
        value_b_list = to_list(value_b)
        return merge_two_lists(value_a_list, value_b_list)
    else:
        # First input node wins
        return value_a


def merge_property_into_existing(
    new_item: dict, existing_item: dict, property_name: str, new_dominates: bool = False
) -> Any:
    dominant_item, secondary_item = (new_item, existing_item) if new_dominates else (existing_item, new_item)
    dominant_value = dominant_item.get(property_name)
    secondary_value = secondary_item.get(property_name)

    # Handle attributes slot specially so we can do nesting at the second level
    if property_name == NODE_ATTRIBUTES or property_name == EDGE_ATTRIBUTES:
        dominant_attributes = dominant_value if dominant_value else dict()
        secondary_attributes = secondary_value if secondary_value else dict()
        source_slots = set(dominant_attributes) | set(secondary_attributes)
        merged_value = {
            source_slot: merge_two_values(dominant_attributes.get(source_slot), secondary_attributes.get(source_slot))
            for source_slot in source_slots
        }
    else:
        merged_value = merge_two_values(dominant_value, secondary_value)

    if property_name == NODE_ID and isinstance(merged_value, list):
        raise ValueError(f"uh oh! ids were merged... shouldn't be possible. {dominant_value}, {secondary_value}")

    existing_item[property_name] = merged_value
