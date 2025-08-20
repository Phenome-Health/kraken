"""
ArangoDB export utilities
Prepares the unified KG for import into ArangoDB
"""

from pathlib import Path
import re
import logging
import json
from ..utils.kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl


def prepare_for_arango_streaming(nodes_file: Path, edges_file: Path, output_dir: Path, config: dict):
    """Prepare unified KG for ArangoDB import using streaming"""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    arango_nodes = output_dir / "arango_nodes.jsonl"
    arango_edges = output_dir / "arango_edges.jsonl"
    
    logging.info("Preparing KG for ArangoDB...")
    
    key_config = config.get('key_field_config', {})
    
    # Process nodes
    node_count = 0
    with open(arango_nodes, 'w') as outfile:
        for node in stream_nodes_from_jsonl(nodes_file):
            # Add ArangoDB _key
            arango_key = create_arango_key(
                node['id'], 
                key_config.get('remove_invalid_chars', True),
                key_config.get('max_length', 254),
                key_config.get('prefix_numbers', True)
            )
            node['_key'] = arango_key
            
            outfile.write(json.dumps(node) + '\n')
            node_count += 1
            
            if node_count % 10000 == 0:
                logging.info(f"Processed {node_count} nodes")
    
    # Process edges
    edge_count = 0
    with open(arango_edges, 'w') as outfile:
        for edge in stream_edges_from_jsonl(edges_file):
            # Add ArangoDB _key
            predicate = edge.get('predicate', 'related_to')
            edge_key = f"{edge['subject']}_{predicate}_{edge['object']}"
            arango_key = create_arango_key(edge_key, True, 254, True)
            edge['_key'] = arango_key
            
            outfile.write(json.dumps(edge) + '\n')
            edge_count += 1
            
            if edge_count % 10000 == 0:
                logging.info(f"Processed {edge_count} edges")
    
    logging.info(f"ArangoDB export complete: {node_count} nodes, {edge_count} edges")
    logging.info(f"Files saved to: {arango_nodes}, {arango_edges}")
    
    return arango_nodes, arango_edges


def create_arango_key(original_key: str, remove_invalid: bool = True,
                      max_length: int = 254, prefix_numbers: bool = True) -> str:
    """Create a valid ArangoDB _key from an original identifier"""
    key = original_key

    if remove_invalid:
        # ArangoDB keys can only contain: a-z, A-Z, 0-9, _, -
        key = re.sub(r'[^a-zA-Z0-9_-]', '_', key)

    if prefix_numbers and key and key[0].isdigit():
        # ArangoDB keys cannot start with a number
        key = f"n_{key}"

    if len(key) > max_length:
        # Truncate but try to keep it unique
        key = key[:max_length - 8] + str(hash(original_key))[-8:]

    return key