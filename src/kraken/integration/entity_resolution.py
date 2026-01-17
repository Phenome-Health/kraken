"""
Entity resolution and graph integration functions
"""

import copy
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import jsonlines

from kraken.utils.constants import (
    AGENT_TYPE,
    CATEGORIES,
    CORE_EDGE_PROPERTIES,
    EQUIVALENT_IDS,
    ID,
    KNOWLEDGE_LEVEL,
    NOT_PROVIDED,
    OBJECT,
    ROOT_CATEGORY,
    SUBJECT,
    SYNONYMS,
)
from kraken.utils.general import create_edge_key, to_list
from kraken.utils.kg_io import (
    get_harmonized_file_paths,
    load_equivalency_mappings,
    remove_file,
    save_to_jsonl,
    stream_edges_from_jsonl,
    stream_nodes_from_jsonl,
)


def integrate_sources(
    source_names: set[str],
    integrated_nodes_path: Path,
    integrated_edges_path: Path,
    harmonized_dir_path: Path,
    config: dict,
):
    """Merge harmonized sources using streaming approach"""
    integrated_dir = integrated_nodes_path.parent
    integrated_dir.mkdir(parents=True, exist_ok=True)

    logging.info("Starting source integration...")

    # Phase 1: Build base equivalency mappings from primary source
    primary_source = config["integration"]["primary_source"]
    primary_nodes_path, _ = get_harmonized_file_paths(primary_source, harmonized_dir_path)
    logging.info(f"Loading equivalency mappings from primary source ({primary_source})")
    equivalency_index = load_equivalency_mappings(primary_nodes_path)
    assert equivalency_index

    # Phase 2: Integrate all nodes, merging as we go
    integrate_nodes(source_names, primary_source, equivalency_index, harmonized_dir_path, integrated_nodes_path, config)

    # Phase 3: Process all edges with node ID resolution (merge edges with the same key -- note aggregator is in key)
    integrate_edges(integrated_edges_path, source_names, equivalency_index, harmonized_dir_path)

    logging.info(f"Integration complete! Unified KG saved to {integrated_dir}")


def integrate_nodes(
    source_names: set[str],
    primary_source: str,
    equivalency_index: dict[str, str],
    harmonized_dir_path: Path,
    integrated_nodes_path: Path,
    config: dict,
):
    integrated_dir = integrated_nodes_path.parent

    # Load the primary source as our starting point
    logging.info(f"Loading {primary_source} nodes as starting point")
    primary_nodes_path, _ = get_harmonized_file_paths(primary_source, harmonized_dir_path)
    current_canonical_nodes = {
        node[ID]: node for node in stream_nodes_from_jsonl(primary_nodes_path)
    }  # canonical_id -> merged_node_data
    assert current_canonical_nodes

    # Figure out what order to integrate sources in (save ones not allowed to merge entities for last)
    sources_config = config["sources"]
    allowed_to_merge_existing = [
        source_name
        for source_name in source_names
        if sources_config[source_name].get("can_merge_existing_nodes") and source_name != primary_source
    ]
    not_allowed_to_merge = [
        source_name
        for source_name in source_names
        if not sources_config[source_name].get("can_merge_existing_nodes") and source_name != primary_source
    ]
    ordered_sources = allowed_to_merge_existing + not_allowed_to_merge
    logging.info(f"Will integrate remaining sources into {primary_source} in this order: {ordered_sources}")

    for source_name in ordered_sources:
        nodes_file, edges_file = get_harmonized_file_paths(source_name, harmonized_dir_path)
        source_allowed_to_merge_nodes = sources_config[source_name].get("can_merge_existing_nodes")

        # Set up logs for non-one-to-one mappings
        one_to_many_log = integrated_dir / f"{source_name}_one_to_many.jsonl"
        one_to_zero_log = integrated_dir / f"{source_name}_one_to_zero.jsonl"
        remove_file(one_to_many_log)
        remove_file(one_to_zero_log)

        logging.info(f"Integrating nodes from {source_name} (can_merge_existing_nodes={source_allowed_to_merge_nodes})")

        for node in stream_nodes_from_jsonl(nodes_file):
            node_id = node[ID]
            node_equiv_ids = node[EQUIVALENT_IDS]

            if source_allowed_to_merge_nodes:
                # Merge all pre-existing canonical nodes referenced by this node's equivalent IDs
                canonical_ids_list = [
                    equivalency_index[equiv_id] for equiv_id in node_equiv_ids if equivalency_index.get(equiv_id)
                ]
                canonical_ids = set(canonical_ids_list)

                if not canonical_ids:
                    # First time seeing this entity in any fashion
                    canonical_id = node[ID]
                    current_canonical_nodes[canonical_id] = node
                    save_to_jsonl([node], one_to_zero_log, mode="a")
                    # Update equivalency index appropriately
                    for equiv_id in node[EQUIVALENT_IDS]:
                        equivalency_index[equiv_id] = canonical_id
                elif len(canonical_ids) == 1:
                    # We have a one-to-one match; merge this node with its canonical node
                    canonical_id = canonical_ids_list[0]
                    existing_canonical_node = current_canonical_nodes[canonical_id]
                    _ = merge_two_nodes(node, existing_canonical_node, equivalency_index, current_canonical_nodes)
                else:
                    # one-to-many match; merge all canonical nodes for this new node into majority canonical node
                    canonical_id_counts = Counter(canonical_ids_list)
                    most_common_canonical_id = canonical_id_counts.most_common(1)[0][0]
                    other_canonical_ids = canonical_ids.difference({most_common_canonical_id})

                    # Log this one-to-many mapping
                    log_item = {
                        "node_id": node_id,
                        "majority_canonical_id": most_common_canonical_id,
                        "other_canonical_ids": list(other_canonical_ids),
                        "node": node,
                    }
                    save_to_jsonl([log_item], one_to_many_log, mode="a")

                    most_common_canonical_node = current_canonical_nodes[most_common_canonical_id]
                    other_canonical_nodes = [current_canonical_nodes[can_id] for can_id in other_canonical_ids]

                    # Merge the new node with the majority canonical node
                    merged_node = merge_two_nodes(
                        node, most_common_canonical_node, equivalency_index, current_canonical_nodes
                    )
                    # Then iteratively merge the other canonical nodes into our merged node
                    for other_existing_canonical_node in other_canonical_nodes:
                        merged_node = merge_two_nodes(
                            other_existing_canonical_node, merged_node, equivalency_index, current_canonical_nodes
                        )
            else:
                # Find the 'majority' canonical ID for this node (not allowed to merge pre-existing canonical nodes)
                canonical_id, new_equiv_ids = find_majority_canonical_id(node, equivalency_index, one_to_many_log)

                if canonical_id in current_canonical_nodes:
                    # Merge with existing canonical node (where new node cannot override existing)
                    # TODO: refine this depending on merging power of pre-existing source(s)?
                    existing_canonical_node = current_canonical_nodes[canonical_id]
                    _ = merge_two_nodes(
                        node,
                        existing_canonical_node,
                        equivalency_index,
                        current_canonical_nodes,
                        new_equiv_ids,
                        new_can_dominate=False,
                    )
                else:
                    # First time seeing this canonical entity
                    current_canonical_nodes[node[ID]] = node
                    save_to_jsonl([node], one_to_zero_log, mode="a")
                    # Update equivalency index appropriately
                    for equiv_id in node[EQUIVALENT_IDS]:
                        equivalency_index[equiv_id] = canonical_id

    logging.info(f"Formed {len(current_canonical_nodes)} merged nodes")

    logging.info("Verifying we have disjoint equivalent_id sets..")
    seen_ids = set()
    for unified_node in current_canonical_nodes.values():
        equiv_ids = set(unified_node[EQUIVALENT_IDS])
        if equiv_ids.intersection(seen_ids):
            logging.error(
                f"Unified node {unified_node[ID]} has equiv IDs present on another unified node(s). "
                f"Overlapping equiv IDs are: {equiv_ids.intersection(seen_ids)}. "
                f"Unified node is: {unified_node}"
            )
            sys.exit(1)
        seen_ids |= equiv_ids

    # Save unified nodes
    save_to_jsonl(current_canonical_nodes.values(), integrated_nodes_path, mode="w")


def integrate_edges(
    integrated_edges_path: Path,
    source_names: set[str],
    equivalency_index: dict[str, str],
    harmonized_dir_path: Path,
):
    integrated_dir = integrated_edges_path.parent
    assert equivalency_index

    all_merged_edges = []
    with jsonlines.open(integrated_edges_path, "w") as writer:
        for source_name in source_names:
            logging.info(f"Processing edges from {source_name}")
            nodes_file, edges_file = get_harmonized_file_paths(source_name, harmonized_dir_path)
            mergers_log = integrated_dir / f"{source_name}_edge_mergers.jsonl"

            # First figure out which edges we're going to need to merge (based on keys)
            edge_key_counts = defaultdict(int)
            for edge in stream_edges_from_jsonl(edges_file):
                # Resolve subject and object to canonical node IDs (needed for accurate keys)
                resolve_to_canonical(edge, equivalency_index)
                key = create_edge_key(edge)
                edge_key_counts[key] += 1
            merged_edges = {key: dict() for key, value in edge_key_counts.items() if value > 1}
            logging.info(f"Identified {len(merged_edges)} {source_name} edges that will be mergers")

            # Then go through and create unified edges
            for edge in stream_edges_from_jsonl(edges_file):
                # Resolve subject and object to canonical IDs
                resolve_to_canonical(edge, equivalency_index)

                # Handle edge merging as necessary
                key = create_edge_key(edge)
                if key in merged_edges:
                    if merged_edges[key]:
                        merge_into_existing_edge(edge, merged_edges[key])  # Add to the merged edge
                    else:
                        merged_edges[key] = edge  # Initiate the merged edge
                else:
                    writer.write(edge)  # No need to merge this edge with others; write it as is

            # Dump all the merged edges for this source to a log for easy review
            logging.info(f"Dumping {len(merged_edges)} merged {source_name} edges to a log..")
            merged_edges_list = list(merged_edges.values())
            save_to_jsonl(merged_edges_list, mergers_log, mode="w")
            all_merged_edges += merged_edges_list

    # Add ALL the edges that had to be merged to our unified edges file (after we closed write mode on the file)
    logging.info(f"Saving {len(all_merged_edges)} total merged edges..")
    save_to_jsonl(all_merged_edges, integrated_edges_path, mode="a")


def find_majority_canonical_id(
    node: dict, equivalency_index: dict[str, str], one_to_many_log: Path
) -> tuple[str, set[str]]:
    """Find canonical ID for this node using equivalency mappings"""
    # Tally up votes for the canonical node from all the equivalent ids
    node_id = node["id"]
    votes = defaultdict(list)
    equiv_ids_without_mappings = set()
    for equiv_id in node[EQUIVALENT_IDS]:
        canonical_id_vote = equivalency_index.get(equiv_id)
        if canonical_id_vote:
            votes[canonical_id_vote].append(equiv_id)
        else:
            equiv_ids_without_mappings.add(equiv_id)
    vote_tallies = {
        canonical_id: len(corresponding_ids)
        + (9 if node_id in corresponding_ids else 0)  # Favor the main node.id (10x the vote)
        for canonical_id, corresponding_ids in votes.items()
    }

    if vote_tallies:
        # Choose the node in the merged graph with the most 'votes' from the equivalent IDs
        canonical_id = max(vote_tallies, key=vote_tallies.get)
        new_equiv_ids = equiv_ids_without_mappings

        # Log if we have a one-to-many mapping
        if len(vote_tallies) > 1:
            log_item = {
                "node_id": node_id,
                "majority_canonical_id": canonical_id,
                "new_equiv_ids": list(new_equiv_ids),
                "vote_tallies": vote_tallies,
                "votes": votes,
                "node": node,
            }
            save_to_jsonl([log_item], one_to_many_log, mode="a")
    else:
        # Can't find a node in the merged graph that this node corresponds to; add it as a new node
        canonical_id = node_id
        new_equiv_ids = set(node[EQUIVALENT_IDS])

    return canonical_id, new_equiv_ids


def merge_two_nodes(
    new_node: dict,
    existing_node: dict,
    equivalency_index: dict[str, str],
    current_canonical_nodes: dict[str, dict],
    new_equiv_ids: set[str] | None = None,
    new_can_dominate: bool = True,
) -> dict[str, Any]:
    """Merge data from new node into existing node (edits in place)"""
    merged_node = copy.deepcopy(existing_node)

    # Figure out whether the new node's values for singular properties should override existing node's
    new_dominates = new_can_dominate and len(new_node[EQUIVALENT_IDS]) > len(existing_node[EQUIVALENT_IDS])

    # Merge any equivalent IDs for this node as appropriate (not necessarily ALL equivalent_ids the source provides,
    #    due to one-to-manys when using majority approach)
    equiv_ids_to_merge = new_equiv_ids if new_equiv_ids is not None else set(new_node[EQUIVALENT_IDS])
    merged_node[EQUIVALENT_IDS] = list(set(existing_node[EQUIVALENT_IDS]) | equiv_ids_to_merge)

    # Only merge in new synonyms if this is a 'full' merge
    if SYNONYMS in new_node and (not new_equiv_ids or len(new_equiv_ids) == len(new_node[EQUIVALENT_IDS])):
        merge_property_into_existing(new_node, merged_node, SYNONYMS)

    # Merge all other properties appropriately
    for property_name, new_value in new_node.items():
        if property_name not in {EQUIVALENT_IDS, SYNONYMS}:  # These are handled specially, above
            merge_property_into_existing(new_node, merged_node, property_name, new_dominates)

    # Remove NamedThing as a category if a more specific category is provided
    if len(merged_node[CATEGORIES]) > 1 and ROOT_CATEGORY in merged_node[CATEGORIES]:
        merged_node[CATEGORIES].remove(ROOT_CATEGORY)

    # Make sure our equivalency index is up to date with any new canonical mappings
    updated_canonical_id = merged_node[ID]
    for equiv_id in merged_node[EQUIVALENT_IDS]:
        equivalency_index[equiv_id] = updated_canonical_id

    # Make sure our canonical nodes map is up to date in light of any changes to canonical ids
    if existing_node[ID] in current_canonical_nodes:
        del current_canonical_nodes[existing_node[ID]]
    if new_node[ID] in current_canonical_nodes:
        del current_canonical_nodes[new_node[ID]]
    current_canonical_nodes[updated_canonical_id] = merged_node

    return merged_node


def merge_into_existing_edge(new_edge: dict, existing_edge: dict):
    # NOTE: If edges are being merged, they must match on all properties included in the edge key

    # Merge knowledge_level, favoring values that aren't not_provided
    if existing_edge[KNOWLEDGE_LEVEL] == NOT_PROVIDED:
        existing_edge[KNOWLEDGE_LEVEL] = new_edge[KNOWLEDGE_LEVEL]

    # Merge agent_type, favoring values that aren't not_provided
    if existing_edge[AGENT_TYPE] == NOT_PROVIDED:
        existing_edge[AGENT_TYPE] = new_edge[AGENT_TYPE]

    # Merge any other properties (all core properties except above 2 are incorporated into key, so must be identical)
    for property_name, value in new_edge.items():
        if property_name not in CORE_EDGE_PROPERTIES:
            merge_property_into_existing(new_edge, existing_edge, property_name)


def resolve_to_canonical(edge: dict, equivalency_index: dict[str, str]):
    subj_id = edge[SUBJECT]
    obj_id = edge[OBJECT]
    if subj_id in equivalency_index and obj_id in equivalency_index:
        edge[SUBJECT] = equivalency_index[edge[SUBJECT]]
        edge[OBJECT] = equivalency_index[edge[OBJECT]]
    else:
        logging.warning(f"Skipping orphan edge: Edge between {subj_id} and {obj_id} is missing equivalency mappings")


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
    if property_name == "attributes":
        dominant_attributes = dominant_value if dominant_value else dict()
        secondary_attributes = secondary_value if secondary_value else dict()
        source_slots = set(dominant_attributes) | set(secondary_attributes)
        merged_value = {
            source_slot: merge_two_values(dominant_attributes.get(source_slot), secondary_attributes.get(source_slot))
            for source_slot in source_slots
        }
    else:
        merged_value = merge_two_values(dominant_value, secondary_value)

    if property_name == ID and isinstance(merged_value, list):
        raise ValueError(f"uh oh! ids were merged... shouldn't be possible. {dominant_value}, {secondary_value}")

    existing_item[property_name] = merged_value
