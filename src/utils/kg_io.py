"""
Knowledge graph I/O utilities
"""
import os
from pathlib import Path
import jsonlines
import logging
from typing import Iterator, Dict, Any


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
