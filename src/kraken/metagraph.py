"""
Metagraph generation utilities for Biolink knowledge graphs
Analyzes node categories, edge predicates, and connectivity patterns
"""

import json
import logging
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from kraken.utils.constants import (
    EDGE_AGENT_TYPE,
    EDGE_KNOWLEDGE_LEVEL,
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_PRIMARY_KS,
    EDGE_SUBJECT,
    EDGE_SUPPORTING_SOURCES,
    NODE_CATEGORIES,
    NODE_EQUIVALENT_IDS,
    NODE_ID,
)
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl


class MetagraphStats:
    """Container for metagraph statistics"""

    def __init__(self, graph_name: str = "unknown", graph_version: str | None = None):
        self.graph_name = graph_name
        self.graph_version = graph_version

        self.node_categories = Counter()  # category -> count
        self.node_prefixes = Counter()  # prefix --> count
        self.total_nodes = 0

        self.total_edges = 0
        self.edge_predicates = Counter()  # predicate -> count
        self.knowledge_sources = Counter()  # primary_knowledge_source + supporting sources -> count
        self.knowledge_levels = Counter()  # knowledge_level -> count
        self.agent_types = Counter()  # agent_type -> count

        self.meta_doubles = Counter()  # (subject_category, object_category) -> count
        self.meta_triples = Counter()  # (subject_cat, predicate, object_cat) -> count

    def to_dict(self) -> dict[str, Any]:
        """Convert stats to dictionary for JSON serialization"""
        return {
            "graph": self.graph_name,
            "version": self.graph_version,
            "summary": {
                "total_nodes": self.total_nodes,
                "total_edges": self.total_edges,
                "unique_node_categories": len(self.node_categories),
                "unique_node_prefixes": len(self.node_prefixes),
                "unique_edge_predicates": len(self.edge_predicates),
                "unique_meta_doubles": len(self.meta_doubles),
                "unique_meta_triples": len(self.meta_triples),
            },
            "node_categories": dict(self.node_categories.most_common()),
            "node_prefixes": dict(self.node_prefixes.most_common()),
            "edge_predicates": dict(self.edge_predicates.most_common()),
            "knowledge_sources": dict(self.knowledge_sources.most_common()),
            "knowledge_levels": dict(self.knowledge_levels.most_common()),
            "agent_types": dict(self.agent_types.most_common()),
            "meta_doubles": {"__".join(double): count for double, count in self.meta_doubles.most_common()},
            "meta_triples": {"__".join(triple): count for triple, count in self.meta_triples.most_common()},
        }


def generate_metagraph_streaming(
    nodes_file: Path, edges_file: Path, graph_name: str, graph_version: str | None
) -> MetagraphStats:
    """Generate metagraph statistics from JSONL files using streaming"""
    logging.info(f"Generating metagraph for {graph_name}")

    stats = MetagraphStats(graph_name, graph_version)

    # Phase 1: Analyze nodes and build category mapping
    categories_map = {}  # node_id -> categories

    logging.info("Analyzing nodes...")
    for node in stream_nodes_from_jsonl(nodes_file):
        node_id = node[NODE_ID]
        categories = node[NODE_CATEGORIES]
        categories_map[node_id] = categories
        for category in categories:
            stats.node_categories[category] += 1

        for equiv_id in node[NODE_EQUIVALENT_IDS]:
            prefix = equiv_id.split(":")[0]
            stats.node_prefixes[prefix] += 1

        stats.total_nodes += 1

    logging.info(f"Found {len(stats.node_categories)} unique node categories")

    if not categories_map:
        logging.error("Categories map is empty.")
        sys.exit(1)

    # Phase 2: Analyze edges
    logging.info("Analyzing edges...")
    for edge in stream_edges_from_jsonl(edges_file):
        subject_id = edge[EDGE_SUBJECT]
        object_id = edge[EDGE_OBJECT]
        predicate = edge[EDGE_PREDICATE]
        if subject_id in categories_map and object_id in categories_map:

            # Collect edge metadata
            if EDGE_PRIMARY_KS in edge:
                stats.knowledge_sources[edge[EDGE_PRIMARY_KS]] += 1
            if EDGE_SUPPORTING_SOURCES in edge:
                for supporting_source in edge[EDGE_SUPPORTING_SOURCES]:
                    stats.knowledge_sources[supporting_source] += 1
            if EDGE_KNOWLEDGE_LEVEL in edge:
                stats.knowledge_levels[edge[EDGE_KNOWLEDGE_LEVEL]] += 1
            if EDGE_AGENT_TYPE in edge:
                stats.agent_types[edge[EDGE_AGENT_TYPE]] += 1

            # Update meta-triple related statistics
            subject_categories = categories_map[subject_id]
            object_categories = categories_map[object_id]
            stats.edge_predicates[predicate] += 1
            for subj_category in subject_categories:
                for obj_category in object_categories:
                    stats.meta_doubles[(subj_category, obj_category)] += 1
                    stats.meta_triples[(subj_category, predicate, obj_category)] += 1

            stats.total_edges += 1
        else:
            logging.warning(f"Orphan edge: Edge between {subject_id} and {object_id} is missing from categories map")

    logging.info(f"Metagraph analysis complete: {stats.total_nodes} nodes, {stats.total_edges} edges")
    return stats


def save_metagraph(stats: MetagraphStats, output_file: Path):
    """Save metagraph statistics to JSON file"""
    logging.info(f"Saving metagraph to {output_file}")

    with open(output_file, "w") as f:
        json.dump(stats.to_dict(), f, indent=2)

    logging.info(f"Metagraph saved: {stats.total_nodes} nodes, {stats.total_edges} edges")


def generate_metagraph_for_source(
    nodes_path: Path, edges_path: Path, output_dir: Path, graph_name: str, graph_version: str = None
) -> Path:
    logging.info(f"Generating metagraph for source {graph_name}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate core statistics
    stats = generate_metagraph_streaming(nodes_path, edges_path, graph_name, graph_version)

    # Save main JSON statistics
    file_name = f"{graph_name}_metagraph_{graph_version}.json" if graph_version else f"{graph_name}_metagraph.json"
    json_file_path = output_dir / file_name
    save_metagraph(stats, json_file_path)

    logging.info(f"Metagraph generated for {graph_name}: {json_file_path}")
    return json_file_path
