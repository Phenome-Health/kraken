"""
Knowledge graph I/O utilities
"""
import os
from pathlib import Path
import json
import jsonlines
import logging
from typing import Iterator, Dict, Any, Optional


def stream_nodes_from_jsonl(nodes_file: Path) -> Iterator[Dict[str, Any]]:
    """Stream nodes from JSONL file without loading into memory"""
    logging.debug(f"Streaming nodes from {nodes_file}")
    
    with jsonlines.open(nodes_file, 'r') as reader:
        for line_num, node in enumerate(reader, 1):
            if line_num % 1000000 == 0:
                logging.info(f"    at {line_num} nodes")
            yield node


def stream_edges_from_jsonl(edges_file: Path) -> Iterator[Dict[str, Any]]:
    """Stream edges from JSONL file without loading into memory"""
    logging.debug(f"Streaming edges from {edges_file}")
    
    with jsonlines.open(edges_file, 'r') as reader:
        for line_num, edge in enumerate(reader, 1):
            if line_num % 5000000 == 0:
                logging.info(f"    at {line_num} edges")
            yield edge


def stream_mixed_jsonl(input_file: Path) -> Iterator[Dict[str, Any]]:
    """Stream items from mixed JSONL file (nodes and edges together)"""
    logging.debug(f"Streaming mixed JSONL from {input_file}")
    
    with jsonlines.open(input_file, 'r') as reader:
        for line_num, item in enumerate(reader, 1):
            if line_num % 1000000 == 0:
                logging.info(f"    at {line_num} items")
            yield item


def save_to_jsonl(items: Iterator[Dict], output_file_path: Path, mode: str = 'w'):
    with jsonlines.open(output_file_path, mode=mode) as writer:
        writer.write_all(items)


def remove_file(file_path: Path):
    if os.path.exists(file_path):
        os.remove(file_path)


def load_node_mappings(nodes_file: Path, key_field: str = 'id') -> Dict[str, Dict]:
    """Load node ID mappings into memory for integration operations"""
    logging.debug(f"Loading node mappings from {nodes_file}")
    
    mappings = {}
    for node in stream_nodes_from_jsonl(nodes_file):
        node_id = node.get(key_field)
        if node_id:
            mappings[node_id] = node
    
    logging.info(f"Loaded {len(mappings)} node mappings")
    return mappings


def load_equivalency_mappings(nodes_file: Path) -> Dict[str, str]:
    """Load equivalency mappings for entity resolution"""
    logging.debug(f"Loading equivalency mappings from {nodes_file}")
    
    equivalencies = {}
    
    for node in stream_nodes_from_jsonl(nodes_file):
        canonical_id = node['id']
        equiv_ids = node['equivalent_ids']

        for equiv_id in equiv_ids:
            equivalencies[equiv_id] = canonical_id
    
    logging.info(f"Loaded equivalencies for {len(equivalencies)} ids")
    return equivalencies


def count_items_in_jsonl(file_path: Path) -> int:
    """Count items in JSONL file without loading into memory"""
    count = 0
    with jsonlines.open(file_path, 'r') as reader:
        for _ in reader:
            count += 1
    return count


def filter_nodes_by_category(nodes: Iterator[Dict], categories: set) -> Iterator[Dict]:
    """Filter nodes by Biolink category"""
    for node in nodes:
        node_categories = node.get('category', [])
        if not isinstance(node_categories, list):
            node_categories = [node_categories]
        
        if any(cat in categories for cat in node_categories):
            yield node


def filter_edges_by_predicate(edges: Iterator[Dict], predicates: set) -> Iterator[Dict]:
    """Filter edges by predicate"""
    for edge in edges:
        if edge.get('predicate') in predicates:
            yield edge


def merge_jsonl_files(input_files: list, output_file: Path):
    """Merge multiple JSONL files into one"""
    logging.info(f"Merging {len(input_files)} files into {output_file}")
    
    total_count = 0
    with jsonlines.open(output_file, 'w') as writer:
        for input_file in input_files:
            if Path(input_file).exists():
                with jsonlines.open(input_file, 'r') as reader:
                    for item in reader:
                        writer.write(item)
                        total_count += 1
    
    logging.info(f"Merged {total_count} items into {output_file}")


# Legacy functions for compatibility (should be gradually replaced)
def save_nodes_edges_as_json(nodes: Iterator[Dict], edges: Iterator[Dict], output_path: Path):
    """Save nodes and edges as traditional JSON format (for compatibility)"""
    logging.debug(f"Saving as JSON to {output_path}")
    
    nodes_list = list(nodes)
    edges_list = list(edges)
    
    data = {
        'nodes': nodes_list,
        'edges': edges_list
    }
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    logging.info(f"Saved {len(nodes_list)} nodes and {len(edges_list)} edges to {output_path}")


def convert_json_to_jsonl(json_file: Path, nodes_output: Path, edges_output: Path):
    """Convert traditional JSON format to JSONL files"""
    logging.info(f"Converting {json_file} to JSONL format")
    
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Save nodes
    if 'nodes' in data:
        save_to_jsonl(iter(data['nodes']), nodes_output)
    
    # Save edges  
    if 'edges' in data:
        save_to_jsonl(iter(data['edges']), edges_output)
