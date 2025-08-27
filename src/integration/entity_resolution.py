"""
Entity resolution and graph integration functions
"""

from collections import defaultdict
from pathlib import Path
import sys
from typing import Dict, Optional, List, Iterator, Set, Tuple, Any
import logging
import jsonlines

from ..utils.kg_io import (
    stream_nodes_from_jsonl, 
    stream_edges_from_jsonl,
    load_equivalency_mappings,
    save_to_jsonl,
    remove_file
)
from ..utils.constants import *
from ..utils.general import create_edge_key


def integrate_sources(harmonized_sources: Dict[str, Dict[str, Path]], output_dir: Path, config: dict):
    """Merge harmonized sources using streaming approach"""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output files
    unified_nodes_path = output_dir / config['unified_output']['nodes']
    unified_edges_path = output_dir / config['unified_output']['edges']
    
    logging.info("Starting streaming source integration...")
    
    # Phase 1: Build equivalency mappings from primary source
    primary_source_name = config.get('primary_source', 'kg2')
    primary_nodes_path = harmonized_sources[primary_source_name]['nodes']
    
    logging.info(f"Loading equivalency mappings from {primary_source_name}")
    equivalency_index = load_equivalency_mappings(primary_nodes_path)
    assert equivalency_index
    logging.info(f"Loading {primary_source_name} nodes as starting point")
    processed_canonical_nodes = {node[ID]: node for node in stream_nodes_from_jsonl(primary_nodes_path)}  # canonical_id -> merged_node_data
    
    # Phase 2: Process all nodes, merging as we go
    for source_name, source_files in harmonized_sources.items():
        if source_name != primary_source_name:  # Already loaded these as the starting point
            # Set up logs for non-one-to-one mappings
            one_to_many_log = output_dir / f"{source_name}_one_to_many.jsonl"
            one_to_zero_log = output_dir / f"{source_name}_one_to_zero.jsonl"
            remove_file(one_to_many_log)
            remove_file(one_to_zero_log)

            logging.info(f"Processing nodes from {source_name}")
            nodes_file = source_files['nodes']
            
            for node in stream_nodes_from_jsonl(nodes_file):
                node_id = node[ID]
                node_equiv_ids = node[EQUIVALENT_IDS]
                
                # Find canonical ID for this node
                canonical_id, new_equiv_ids = find_canonical_id(node_id, node_equiv_ids, equivalency_index, one_to_many_log, node)
                
                if canonical_id in processed_canonical_nodes:
                    # Merge with existing canonical node
                    existing_canonical_node = processed_canonical_nodes[canonical_id]
                    merge_into_existing_node(node, existing_canonical_node, new_equiv_ids)
                else:
                    # First time seeing this canonical entity
                    processed_canonical_nodes[node[ID]] = node
                    save_to_jsonl([node], one_to_zero_log, mode='a')
                
                # Update our equivalency index with any new canonical mappings
                for equiv_id in processed_canonical_nodes[canonical_id][EQUIVALENT_IDS]:
                    equivalency_index[equiv_id] = canonical_id

    
    logging.info(f"Formed {len(processed_canonical_nodes)} merged nodes")

    logging.info(f"Verifying we have disjoint equivalent_id sets..")
    seen_ids = set()
    for unified_node in processed_canonical_nodes.values():
        equiv_ids = set(unified_node[EQUIVALENT_IDS])
        if equiv_ids.intersection(seen_ids):
            logging.error(f"Unified node {unified_node[ID]} has equiv IDs present on another unified node(s). "
                          f"Overlapping equiv IDs are: {equiv_ids.intersection(seen_ids)}")
            sys.exit(1)
        seen_ids |= equiv_ids

    # Save unified nodes
    save_to_jsonl(processed_canonical_nodes.values(), unified_nodes_path, mode='w')
    
    # Phase 3: Process all edges with node ID resolution (merge edges with the same key -- note aggregator is in key)

    with jsonlines.open(unified_edges_path, 'w') as writer:
        for source_name, source_files in harmonized_sources.items():
            logging.info(f"Processing edges from {source_name}")
            edges_file = source_files['edges']

            # First figure out which edges we're going to need to merge (based on keys)
            edge_key_counts = defaultdict(int)
            for edge in stream_edges_from_jsonl(edges_file):
                # Resolve subject and object to canonical node IDs (needed for accurate keys)
                resolve_to_canonical(edge, equivalency_index)
                key = create_edge_key(edge)
                edge_key_counts[key] += 1
            assert edge_key_counts
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

        # Now dump all the edges from this source that had to be merged
        logging.info(f"Saving {len(merged_edges)} merged {source_name} edges..")
        save_to_jsonl(merged_edges.values(), unified_edges_path, mode='a')
    
    logging.info(f"Integration complete! Unified KG saved to {output_dir}")


def find_canonical_id(node_id: str, equiv_ids: List[str], equivalency_index: Dict[str, str], one_to_many_log: Path, node: dict) -> Tuple[str, Set[str]]:
    """Find canonical ID for this node using equivalency mappings"""
    # Tally up votes for the canonical node from all the equivalent ids
    votes = defaultdict(list)
    no_mappings = set()
    for equiv_id in equiv_ids:
        canonical_id = equivalency_index.get(equiv_id)
        if canonical_id:
            votes[canonical_id].append(equiv_id)
        else:
            no_mappings.add(equiv_id)
    vote_tallies = {canonical_id: len(corresponding_ids) + (9 if node_id in corresponding_ids else 0)  # Favor the main node.id (10x the vote)
                    for canonical_id, corresponding_ids in votes.items()}
    
    if vote_tallies:
        # Choose the node in the merged graph with the most 'votes' from the equivalent IDs
        canonical_id = max(vote_tallies, key=vote_tallies.get)
        new_equiv_ids = no_mappings

        # Log if we have a one-to-many mapping
        if len(vote_tallies) > 1:
            log_item = {'node_id': node_id, 'chosen_canonical_id': canonical_id, 'equivalent_ids': equiv_ids,
                        'new_equiv_ids': list(new_equiv_ids), 'votes': votes, 'vote_tallies': vote_tallies, 'node': node}
            save_to_jsonl([log_item], one_to_many_log, mode='a')
    else:
        # Can't find a node in the merged graph that this node corresponds to; add it as a new node
        canonical_id = node_id
        new_equiv_ids = equiv_ids
    
    return canonical_id, new_equiv_ids


def merge_into_existing_node(new_node: dict, existing_node: dict, new_equiv_ids: Set[str]):
    """Merge data from new node into existing node (edits in place)"""

    # Add any 'new' equivalent IDs for this node (not necessarily ALL equivalent_ids the source provides, due to one-to-manys)
    existing_node[EQUIVALENT_IDS] = list(set(existing_node[EQUIVALENT_IDS]) | new_equiv_ids)
    del new_node[EQUIVALENT_IDS]  # We don't want any one-to-manys that lost the vote appearing here

    # Merge other list properties
    existing_node[CATEGORIES] = merge_two_list_properties(existing_node, new_node, CATEGORIES)
    existing_node[PROVIDED_BY] =  merge_two_list_properties(existing_node, new_node, PROVIDED_BY)
    if SYNONYMS in existing_node or SYNONYMS in new_node:
        existing_node[SYNONYMS] = merge_two_list_properties(existing_node, new_node, SYNONYMS)
    if new_node.get(NAME):  # Add the new node's name as a synonym for the merged node
        existing_node[SYNONYMS] = list(set(existing_node.get(SYNONYMS, [])) | {new_node[NAME]})
    
    # Remove NamedThing as a category if a more specific category is provided
    if len(existing_node[CATEGORIES]) > 1 and ROOT_CATEGORY in existing_node[CATEGORIES]:
        existing_node[CATEGORIES].remove(ROOT_CATEGORY)
    
    # Merge any other properties appropriately
    for property_name, value in new_node.items():
        if property_name not in CORE_NODE_PROPERTIES:
            if isinstance(value, list):
                if any(isinstance(item, dict) for item in value):
                    existing_node[property_name] = concatenate_two_list_properties(existing_node, new_node, property_name)
                else:
                    existing_node[property_name] = merge_two_list_properties(existing_node, new_node, property_name)
            else:
                # First come first serve
                if existing_node.get(property_name) is None:
                    existing_node[property_name] = value


def merge_into_existing_edge(new_edge: dict, existing_edge: dict):
    # NOTE: If edges are being merged, they must match on all properties included in the edge key

    # Merge knowledge_level
    if new_edge.get(KNOWLEDGE_LEVEL):
        if not existing_edge.get(KNOWLEDGE_LEVEL) or existing_edge[KNOWLEDGE_LEVEL] == UNKNOWN_KNOWLEDGE_LEVEL:
            existing_edge[KNOWLEDGE_LEVEL] = new_edge[KNOWLEDGE_LEVEL]

    # Merge agent_type
    if new_edge.get(AGENT_TYPE):
        if not existing_edge.get(AGENT_TYPE) or existing_edge[AGENT_TYPE] == UNKNOWN_AGENT_TYPE:
            existing_edge[AGENT_TYPE] = new_edge[AGENT_TYPE]

    # Merge any other properties (all core properties except above 2 are incorporated into key, so must be identical)
    for property_name, value in new_edge.items():
        if property_name not in CORE_EDGE_PROPERTIES:
            if isinstance(value, list):
                if any(isinstance(item, dict) for item in value):
                    existing_edge[property_name] = concatenate_two_list_properties(existing_edge, new_edge, property_name)
                else:
                    existing_edge[property_name] = merge_two_list_properties(existing_edge, new_edge, property_name)
            else:
                # First come first serve
                if existing_edge.get(property_name) is None:
                    existing_edge[property_name] = value


def merge_two_list_properties(node_a: dict, node_b: dict, property_name: str) -> List[Any]:
    # Merges two list properties, retaining distinct values
    return list(set(node_a.get(property_name, [])) | set(node_b.get(property_name, [])))


def concatenate_two_list_properties(node_a: dict, node_b: dict, property_name: str) -> List[Any]:
    # Concatenates two list properties (does not check for uniqueness of items)
    return node_a.get(property_name, []) + node_b.get(property_name, [])


def resolve_to_canonical(edge: dict, equivalency_index: Dict[str, str]):
    subj_id = edge[SUBJECT]
    obj_id = edge[OBJECT]
    if subj_id in equivalency_index and obj_id in equivalency_index:
        edge[SUBJECT] = equivalency_index[edge[SUBJECT]]
        edge[OBJECT] = equivalency_index[edge[OBJECT]]
    else:
        logging.warning(f"Skipping orphan edge: Edge between {subj_id} and {obj_id} is missing equivalency mappings")
