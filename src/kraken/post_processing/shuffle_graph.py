"""
Create a 'shuffled' version of KRAKEN by permuting node IDs in edges.

All node IDs are randomly shuffled via a permutation map (by category), then each edge's
subject and object are replaced with their shuffled counterparts. The shuffled
edges are written to a new file alongside the original.

Usage:
    python shuffle_kraken.py --nodes <nodes.jsonl> --edges <edges.jsonl> [--seed 42]
"""

import argparse
import logging
import random
import time
from collections import defaultdict
from pathlib import Path

import jsonlines

from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl
from kraken.utils.logging_config import setup_logging


def create_permutation_map(nodes_path: Path, rng: random.Random) -> dict[str, str]:
    """Build a permutation map from node IDs to shuffled node IDs, shuffling within categories."""
    logging.info(f"Loading node IDs by category from {nodes_path}...")

    category_to_nodes: dict[str, list[str]] = defaultdict(list)
    total_nodes = 0

    for node in stream_nodes_from_jsonl(nodes_path):
        categories = node["categories"]
        category = rng.choice(categories)
        category_to_nodes[category].append(node["id"])
        total_nodes += 1

    logging.info(f"Loaded {total_nodes:,} nodes across {len(category_to_nodes):,} categories")

    perm_map: dict[str, str] = {}
    for category, nodes in category_to_nodes.items():
        indices = list(range(len(nodes)))
        rng.shuffle(indices)
        for i, idx in enumerate(indices):
            perm_map[nodes[i]] = nodes[idx]

    logging.info(f"Built permutation map ({len(perm_map):,} entries)")
    return perm_map


def shuffle_edges(edges_path: Path, perm_map: dict[str, str], output_dir: Path) -> Path:
    """Write shuffled edges to a new file, replacing subject/object via the permutation map."""
    # TODO: make work with passing in outpudir from orchestrator?
    output_file_name = edges_path.stem + "_shuffled" + edges_path.suffix
    output_path = output_dir / output_file_name
    logging.info(f"Writing shuffled edges to {output_path}...")

    edge_count = 0
    with jsonlines.open(output_path, "w") as writer:
        for edge in stream_edges_from_jsonl(edges_path):
            subj = edge["subject"]
            obj = edge["object"]
            edge["subject"] = perm_map[subj]
            edge["object"] = perm_map[obj]
            writer.write(edge)
            edge_count += 1

    logging.info(f"Wrote {edge_count:,} shuffled edges total")
    return output_path


def shuffle_graph(nodes_path: Path, edges_path: Path, output_dir: Path, seed: int = 42):
    start = time.time()
    logging.info(f"Shuffling graph with seed={seed}")

    rng = random.Random(seed)
    perm_map = create_permutation_map(nodes_path, rng)
    output_path = shuffle_edges(edges_path, perm_map, output_dir)

    elapsed = time.time() - start
    logging.info(f"Done! Wrote shuffled edges to {output_path}")
    logging.info(f"Total time: {elapsed:.1f}s ({elapsed / 60:.1f} min)")


def main():
    setup_logging()

    parser = argparse.ArgumentParser(description="Create a shuffled version of KRAKEN")
    parser.add_argument("--nodes", type=Path, required=True, help="Path to nodes JSONL file")
    parser.add_argument("--edges", type=Path, required=True, help="Path to edges JSONL file")
    parser.add_argument("--output", type=Path, required=True, help="Path to output directory")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    args = parser.parse_args()

    shuffle_graph(args.nodes, args.edges, args.output, args.seed)


if __name__ == "__main__":
    main()
