"""
Entity resolution and graph integration functions
"""

from collections import defaultdict
from pathlib import Path
import sys
from typing import Dict, Optional, List, Iterator, Set, Tuple
import logging
import jsonlines

from ..utils.kg_io import (
    stream_nodes_from_jsonl, 
    stream_edges_from_jsonl,
    load_equivalency_mappings,
    save_to_jsonl,
    remove_file
)

LIST_PROPERTIES = ["synonyms", "provided_by", "categories"]


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
    processed_canonical_nodes = {node['id']: node for node in stream_nodes_from_jsonl(primary_nodes_path)}  # canonical_id -> merged_node_data
    
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
                node_id = node['id']
                node_equiv_ids = node['equivalent_ids']
                
                # Find canonical ID for this node
                canonical_id, new_equiv_ids = find_canonical_id(node_id, node_equiv_ids, equivalency_index, one_to_many_log, node)
                
                if canonical_id in processed_canonical_nodes:
                    # Merge with existing canonical node
                    existing_canonical_node = processed_canonical_nodes[canonical_id]
                    merge_into_existing_node(node, existing_canonical_node, new_equiv_ids, source_name)
                else:
                    # First time seeing this canonical entity
                    processed_canonical_nodes[node['id']] = node
                    save_to_jsonl([node], one_to_zero_log, mode='a')
                
                # Update our equivalency index with any new canonical mappings
                for equiv_id in processed_canonical_nodes[canonical_id]['equivalent_ids']:
                    equivalency_index[equiv_id] = canonical_id

    
    logging.info(f"Formed {len(processed_canonical_nodes)} merged nodes")

    logging.info(f"Verifying we have disjoint equivalent_id sets..")
    seen_ids = set()
    for unified_node in processed_canonical_nodes.values():
        equiv_ids = set(unified_node['equivalent_ids'])
        if equiv_ids.intersection(seen_ids):
            logging.error(f"Unified node {unified_node['id']} has equiv IDs present on another unified node(s). "
                          f"Overlapping equiv IDs are: {equiv_ids.intersection(seen_ids)}")
            sys.exit(1)
        seen_ids |= equiv_ids

    # Save unified nodes
    save_to_jsonl(processed_canonical_nodes.values(), unified_nodes_path)
    
    # Phase 3: Process all edges with node ID resolution (no merging)

    with jsonlines.open(unified_edges_path, 'w') as writer:
        for source_name, source_files in harmonized_sources.items():
            logging.info(f"Processing edges from {source_name}")
            edges_file = source_files['edges']
            
            for edge in stream_edges_from_jsonl(edges_file):
                # Resolve subject and object to canonical IDs
                subj_id = edge['subject']
                obj_id = edge['object']
                if subj_id in equivalency_index and obj_id in equivalency_index:
                    edge['subject'] = equivalency_index[edge['subject']]
                    edge['object'] = equivalency_index[edge['object']]
                else:
                    logging.warning(f"Skipping oprhan edge: Edge between {subj_id} and {obj_id} is missing equivalency mappings")
                writer.write(edge)
    
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


def merge_into_existing_node(new_node: dict, existing_node: dict, new_equiv_ids: Set[str], source_name: str):
    """Merge data from new node into existing node (edits in place)"""
    # Add any new equivalent IDs for this node (not necessarily all equivalent_ids the source provides)
    existing_node['equivalent_ids'] = list(set(existing_node['equivalent_ids']) | new_equiv_ids)
    del new_node['equivalent_ids']  # We don't want any one-to-manys that lost the vote appearing here

    # Merge other list properties
    for property_name in LIST_PROPERTIES:
        if property_name in existing_node or property_name in new_node:
            existing_node[property_name] = list(set(existing_node.get(property_name, [])) | set(new_node.get(property_name, [])))
    
    # Remove NamedThing as a category if a more specific category is provided
    if len(existing_node['categories']) > 1 and 'biolink:NamedThing' in existing_node['categories']:
        existing_node['categories'].remove('biolink:NamedThing')
    
    # Merge other properties (simple first-wins strategy for now)
    for key, value in new_node.items():
        if key == 'id':
            # Keep existing canonical ID
            continue
        elif key not in existing_node or existing_node[key] is None:
            # Add new property
            existing_node[key] = value
        # For conflicts, keep existing value (could be made more sophisticated)
