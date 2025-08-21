"""
Entity resolution and graph integration functions
"""

from pathlib import Path
from typing import Dict, Optional, List, Iterator
import logging
import jsonlines

from ..utils.kg_io import (
    stream_nodes_from_jsonl, 
    stream_edges_from_jsonl,
    load_equivalency_mappings,
    save_nodes_to_jsonl
)
from ..utils.metagraph import generate_metagraph_for_source, compare_metagraphs


def integrate_sources(harmonized_sources: Dict[str, Dict[str, Path]], output_dir: Path, config: dict):
    """Merge harmonized sources using streaming approach"""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output files
    unified_nodes_path = output_dir / "unified_nodes.jsonl"
    unified_edges_path = output_dir / "unified_edges.jsonl"
    
    logging.info("Starting streaming source integration...")
    
    # Phase 1: Build equivalency mappings from primary source
    primary_source_name = config.get('primary_source', 'rtx-kg2')
    primary_nodes_path = harmonized_sources[primary_source_name]['nodes']
    
    logging.info(f"Loading equivalency mappings from {primary_source_name}")
    equivalency_index = load_equivalency_mappings(primary_nodes_path)
    processed_canonical_nodes = {node['id']: node for node in stream_nodes_from_jsonl(primary_nodes_path)}  # canonical_id -> merged_node_data
    
    # Phase 2: Process all nodes, merging as we go

    for source_name, source_files in harmonized_sources.items():
        if source_name != primary_source_name:  # Already loaded these as the starting point
            logging.info(f"Processing nodes from {source_name}")
            nodes_file = source_files['nodes']
            
            for node in stream_nodes_from_jsonl(nodes_file):
                node_id = node['id']
                node_equiv_ids = node['equivalent_ids']
                
                # Find canonical ID for this node
                canonical_id = find_canonical_id(node_id, node_equiv_ids, equivalency_index)
                
                if canonical_id in processed_canonical_nodes:
                    # Merge with existing canonical node
                    existing_canonical_node = processed_canonical_nodes[canonical_id]
                    merge_into_existing_node(node, existing_canonical_node)
                else:
                    # First time seeing this canonical entity
                    processed_canonical_nodes[canonical_id] = node
                    processed_canonical_nodes[canonical_id]['id'] = canonical_id
                    # Update our equivalency index with these new canonical mappings
                    for equiv_id in node_equiv_ids:
                        equivalency_index[equiv_id] = canonical_id

    
    logging.info(f"Formed {len(processed_canonical_nodes)} merged nodes")
    
    # Save unified nodes
    save_nodes_to_jsonl(processed_canonical_nodes.values(), unified_nodes_path)
    
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
    
    # Generate metagraph for unified result
    if config.get('generate_metagraph', True):
        # Store unified metagraphs in artifacts/metagraphs/unified/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "unified"
        
        metagraph_config = config.get('metagraph_config', {
            'generate_summaries': True,
            'generate_cytoscape': True,
            'generate_html_viewer': True,
            'cytoscape_thresholds': [1, 5, 10, 25]  # Additional threshold for unified graph
        })
        
        unified_metagraph_files = generate_metagraph_for_source(
            unified_nodes_path, unified_edges_path, metagraph_dir, "unified", metagraph_config
        )
        logging.info("Unified metagraph generated")
        
        # Compare with source metagraphs if they exist
        source_metagraphs = []
        for source_name in harmonized_sources.keys():
            source_metagraph = artifacts_root / "metagraphs" / "harmonized" / source_name / f"{source_name}_metagraph.json"
            if source_metagraph.exists():
                source_metagraphs.append(source_metagraph)
        
        if source_metagraphs:
            # Find the main JSON file from unified metagraph
            unified_json = next((f for f in unified_metagraph_files if f.name.endswith('_metagraph.json')), None)
            if unified_json:
                source_metagraphs.append(unified_json)
                comparison_file = metagraph_dir / "metagraph_comparison.json"
                compare_metagraphs(source_metagraphs, comparison_file)
                logging.info("Metagraph comparison generated")
    
    return unified_nodes_path, unified_edges_path


def find_canonical_id(node_id: str, equiv_ids: List[str], equivalency_index: Dict[str, set]) -> str:
    """Find canonical ID for this node using equivalency mappings"""
    # Check if this node or any equivalent ID is in our index
    all_ids = [node_id] + equiv_ids
    
    for id_val in all_ids:
        if id_val in equivalency_index:
            # Return the first (canonical) ID from the equivalency set
            return equivalency_index[id_val]
    
    # Not found in index, use the original ID as canonical
    return node_id


def merge_into_existing_node(new_node: dict, existing_node: dict):
    """Merge data from new node into existing node (edits in place)"""
    # Merge equivalent_ids
    existing_equivs = set(existing_node.get('equivalent_ids', []))
    new_equivs = set(new_node.get('equivalent_ids', []))
    existing_node['equivalent_ids'] = list(existing_equivs | new_equivs)

    # Merge sources
    existing_node['provided_by'] = existing_node.get('provided_by', []) + new_node.get('provided_by', [])

    # Merge other properties (simple first-wins strategy for now)
    for key, value in new_node.items():
        if key == 'id':
            # Keep existing canonical ID
            continue
        elif key not in existing_node or existing_node[key] is None:
            # Add new property
            existing_node[key] = value
        # For conflicts, keep existing value (could be made more sophisticated)
