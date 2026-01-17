import logging
import random
from pathlib import Path

import jsonlines

from ..utils.general import create_edge_key
from ..utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl, save_to_jsonl


def create_test_kg_files(nodes_path: Path, edges_path: Path, output_dir: Path, num_edges: int = 1000):
    logging.info(f"Creating test versions of the KG JSON Lines files at {nodes_path} and {edges_path}...")
    logging.info(f"Number of edges to includes is {num_edges}")

    logging.info(f"Loading all edge IDs")
    edge_ids = [create_edge_key(edge) for edge in stream_edges_from_jsonl(edges_path)]
    logging.info(f"Loaded {len(edge_ids)} edge IDs")

    logging.info(f"Randomly selecting {num_edges} of those edges")
    test_edge_ids = set(random.sample(edge_ids, num_edges))
    test_edges = [edge for edge in stream_edges_from_jsonl(edges_path) if create_edge_key(edge) in test_edge_ids]
    logging.info(f"Grabbed {len(test_edges)} random test edges")
    assert len(test_edges) == len(test_edge_ids)

    logging.info(f"Grabbing nodes that are used by those edges")
    node_ids_used = {edge['subject'] for edge in test_edges}.union({edge['object'] for edge in test_edges})
    test_nodes = [node for node in stream_nodes_from_jsonl(nodes_path) if node['id'] in node_ids_used]
    logging.info(f"Grabbed the {len(test_nodes)} nodes used by the selected edges.")
    assert len(node_ids_used) == len(test_nodes)

    test_nodes_file_name = nodes_path.name.replace('.jsonl', '_test.jsonl')
    test_edges_file_name = edges_path.name.replace('.jsonl', '_test.jsonl')
    test_nodes_path = output_dir / test_nodes_file_name
    test_edges_path = output_dir / test_edges_file_name
    logging.info(f"Saving {len(test_nodes)} test nodes to {test_nodes_path}")
    logging.info(f"Saving {len(test_edges)} test edges to {test_edges_path}")
    save_to_jsonl(test_nodes, test_nodes_path, mode='w')
    save_to_jsonl(test_edges, test_edges_path, mode='w')
