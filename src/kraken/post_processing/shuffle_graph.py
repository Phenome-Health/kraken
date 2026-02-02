"""
Shuffle graph edges for negative control generation.

This script creates a 'shuffled' version of the KG where subject/object IDs
are randomly permuted within the same Biolink category, and ensures no shuffled
edge already exists in the original graph.

Usage:
    python shuffle_graph.py --nodes nodes.jsonl --edges edges.jsonl --output-dir shuffled/
"""

import argparse
import logging
import random
from collections import defaultdict
from itertools import chain
from pathlib import Path

import jsonlines

from kraken.utils.general import create_edge_key
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl
from kraken.utils.logging_config import setup_logging

setup_logging()


def build_category_maps(
    nodes_path: Path,
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    """
    Build mappings between node IDs and their Biolink categories.

    Returns:
        category_to_node_list: Maps each category to a list of node IDs (for efficient sampling)
        node_id_to_categories: Maps each node ID to all its categories
    """
    category_to_node_ids: dict[str, set[str]] = defaultdict(set)
    node_id_to_categories: dict[str, list[str]] = {}

    for node in stream_nodes_from_jsonl(nodes_path):
        node_id = node["id"]
        categories = node.get("categories", [])
        if not categories:
            # Fallback if no categories - shouldn't happen in well-formed data
            categories = ["biolink:NamedThing"]

        node_id_to_categories[node_id] = categories

        # Register this node under ALL its categories
        for category in categories:
            category_to_node_ids[category].add(node_id)

    # Convert to lists for efficient random sampling
    category_to_node_list: dict[str, list[str]] = {cat: list(nodes) for cat, nodes in category_to_node_ids.items()}

    return category_to_node_list, node_id_to_categories


def collect_existing_edge_keys(edges_path: Path) -> tuple[set[str], int]:
    """
    Collect all edge keys from the original graph.
    """
    existing_keys: set[str] = set()
    edge_count = 0

    for edge in stream_edges_from_jsonl(edges_path):
        edge_key = create_edge_key(edge)
        existing_keys.add(edge_key)
        edge_count += 1

    return existing_keys, edge_count


def get_shuffled_id(
    original_id: str,
    node_id_to_categories: dict[str, list[str]],
    category_to_node_list: dict[str, list[str]],
    rng: random.Random,
) -> str:
    """
    Get a random node ID that shares at least one category with the original.
    """
    categories = node_id_to_categories.get(original_id)
    if categories is None:
        raise ValueError(f"Node {original_id} not present in node_id_to_categories map")

    # Shuffle categories so we try them in random order
    rng.shuffle(categories)

    for category in categories:
        candidates = category_to_node_list[category]

        if len(candidates) < 2:
            continue

        # Try a few times to get something other than the original
        for _ in range(10):
            chosen = rng.choice(candidates)
            if chosen != original_id:
                return chosen

    # No category had valid candidates
    return original_id


def shuffle_edge(
    edge: dict,
    node_id_to_categories: dict[str, list[str]],
    category_to_node_list: dict[str, list[str]],
    existing_keys: set[str],
    shuffled_keys: set[str],
    rng: random.Random,
    max_attempts: int = 100,
) -> tuple[str, dict]:
    """
    Shuffle an edge's subject and object IDs in place.

    Shuffles subject and object IDs (within overlapping categories) and ensures:
    - The resulting edge doesn't already exist in original or shuffled graph
    - No self-loops (subject != object)

    Returns:
        Shuffled edge (dict) and its key (string)
    """
    original_subject = edge["subject"]
    original_object = edge["object"]

    for _ in range(max_attempts):
        # Shuffle both subject and object
        new_subject = get_shuffled_id(original_subject, node_id_to_categories, category_to_node_list, rng)
        new_object = get_shuffled_id(original_object, node_id_to_categories, category_to_node_list, rng)

        # Skip if we got the exact same edge back
        if new_subject == original_subject and new_object == original_object:
            continue

        # Skip self-loops
        if new_subject == new_object:
            continue

        # Temporarily update edge to compute new shuffled key
        edge["subject"] = new_subject
        edge["object"] = new_object

        shuffled_key = create_edge_key(edge)

        if shuffled_key not in existing_keys and shuffled_key not in shuffled_keys:
            return shuffled_key, edge

    # Failed - restore original values and return empty
    edge["subject"] = original_subject
    edge["object"] = original_object
    return "", dict()


def shuffle_graph(
    nodes_path: Path,
    edges_path: Path,
    output_dir: Path,
    seed: int = 42,
) -> None:
    """
    Main function to create shuffled graph.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_edges_path = output_dir / edges_path.name.replace(".jsonl", "_shuffled.jsonl")

    # Initialize random number generator
    rng = random.Random(seed)

    # Step 1: Build category maps from nodes
    logging.info("Building category maps from nodes...")
    category_to_node_list, node_id_to_categories = build_category_maps(nodes_path)
    logging.info(f"  Found {len(node_id_to_categories)} nodes across {len(category_to_node_list)} categories")

    # Step 2: Collect existing edge keys
    logging.info("Collecting existing edge keys...")
    existing_keys, existing_edge_count = collect_existing_edge_keys(edges_path)
    logging.info(f"  Found {len(existing_keys)} existing edge keys, from {existing_edge_count} existing edges")

    # Step 3: Shuffle edges and write to output
    logging.info("Shuffling edges...")
    shuffled_keys: set[str] = set()
    edges_written = 0
    edges_skipped = 0

    with jsonlines.open(output_edges_path, "w") as writer:
        for edge in stream_edges_from_jsonl(edges_path):
            shuffled_key, shuffled_edge = shuffle_edge(
                edge,
                node_id_to_categories,
                category_to_node_list,
                existing_keys,
                shuffled_keys,
                rng,
            )

            if shuffled_edge:
                # Track this edge to avoid duplicates in shuffled output
                shuffled_keys.add(shuffled_key)

                # Write to shuffled graph file
                writer.write(shuffled_edge)
                edges_written += 1
            else:
                edges_skipped += 1

            # Progress indicator
            total_processed = edges_written + edges_skipped
            if total_processed % 10_000_000 == 0:
                logging.info(
                    f"  Processed {total_processed:,} edges ({edges_written} written, {edges_skipped} skipped)..."
                )

    logging.info(f"Done! Wrote {edges_written:,} shuffled edges to {output_edges_path}")
    if edges_skipped > 0:
        logging.warning(f"  Skipped {edges_skipped:,} edges (couldn't find valid shuffle)")

    assert (edges_written + edges_skipped) == existing_edge_count


def main():
    parser = argparse.ArgumentParser(description="Create shuffled version of KG for negative control generation")
    parser.add_argument(
        "--nodes",
        type=Path,
        required=True,
        help="Path to nodes.jsonl file",
    )
    parser.add_argument(
        "--edges",
        type=Path,
        required=True,
        help="Path to edges.jsonl file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for shuffled files",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )

    args = parser.parse_args()

    shuffle_graph(
        nodes_path=args.nodes,
        edges_path=args.edges,
        output_dir=args.output_dir,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
