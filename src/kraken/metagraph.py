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
    EDGE_AGGREGATOR_KS,
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


def _branch_total(node: Any) -> int:
    """Total of all leaf counts under a nested node (or the count itself if a leaf)."""
    return sum(_branch_total(child) for child in node.values()) if isinstance(node, dict) else node


def _sort_branch(node: dict) -> dict:
    """Recursively order each level of a nested dict by descending subtree total,
    so the heaviest branches (and, at the leaves, the most common counts) come first."""
    return {
        key: (_sort_branch(child) if isinstance(child, dict) else child)
        for key, child in sorted(node.items(), key=lambda item: _branch_total(item[1]), reverse=True)
    }


def _nest_and_sort_counts(counter: Counter) -> dict:
    """Nest a tuple-keyed Counter (e.g. (subject, predicate, object) -> count) into nested
    dicts keyed by each tuple element, ordered by descending count at every level."""
    nested: dict = {}
    for key_tuple, count in counter.items():
        level = nested
        for part in key_tuple[:-1]:
            level = level.setdefault(part, {})
        level[key_tuple[-1]] = count
    return _sort_branch(nested)


class MetagraphStats:
    """Container for metagraph statistics"""

    def __init__(
        self,
        graph_name: str = "unknown",
        graph_version: str | None = None,
        source_versions: dict[str, str | None] | None = None,
        biolink_version: str | None = None,
    ):
        self.graph_name = graph_name
        self.graph_version = graph_version
        self.biolink_version = biolink_version  # Biolink Model version the graph conforms to
        self.source_versions = source_versions  # {source_name: version} of ingested sources

        self.node_categories = Counter()  # category -> count
        self.node_prefixes = Counter()  # prefix --> count
        self.total_nodes = 0

        self.total_edges = 0
        self.edge_predicates = Counter()  # predicate -> count
        # Provenance is tracked per Biolink role, not conflated into one bucket:
        self.primary_knowledge_sources = Counter()  # primary_knowledge_source -> count
        self.aggregator_knowledge_sources = Counter()  # aggregator_knowledge_source -> count
        self.supporting_data_sources = Counter()  # supporting_data_sources -> count
        self.knowledge_levels = Counter()  # knowledge_level -> count
        self.agent_types = Counter()  # agent_type -> count
        self.klat_joint = Counter()  # (knowledge_level, agent_type) -> count

        self.meta_doubles = Counter()  # (subject_category, object_category) -> count
        self.meta_triples = Counter()  # (subject_cat, predicate, object_cat) -> count

    def to_dict(self) -> dict[str, Any]:
        """Convert stats to dictionary for JSON serialization"""
        return {
            "graph": self.graph_name,
            "version": self.graph_version,
            "biolink_version": self.biolink_version,
            "source_versions": self.source_versions,
            "summary": {
                "total_nodes": self.total_nodes,
                "total_edges": self.total_edges,
                "unique_node_categories": len(self.node_categories),
                "unique_node_prefixes": len(self.node_prefixes),
                "unique_edge_predicates": len(self.edge_predicates),
                "unique_primary_knowledge_sources": len(self.primary_knowledge_sources),
                "unique_aggregator_knowledge_sources": len(self.aggregator_knowledge_sources),
                "unique_supporting_data_sources": len(self.supporting_data_sources),
                "unique_meta_doubles": len(self.meta_doubles),
                "unique_meta_triples": len(self.meta_triples),
                "unique_klat_combinations": len(self.klat_joint),
            },
            "node_categories": dict(self.node_categories.most_common()),
            "node_prefixes": dict(self.node_prefixes.most_common()),
            "edge_predicates": dict(self.edge_predicates.most_common()),
            "primary_knowledge_sources": dict(self.primary_knowledge_sources.most_common()),
            "aggregator_knowledge_sources": dict(self.aggregator_knowledge_sources.most_common()),
            "supporting_data_sources": dict(self.supporting_data_sources.most_common()),
            "knowledge_levels": dict(self.knowledge_levels.most_common()),
            "agent_types": dict(self.agent_types.most_common()),
            "klat_joint": _nest_and_sort_counts(self.klat_joint),
            "meta_doubles": _nest_and_sort_counts(self.meta_doubles),
            "meta_triples": _nest_and_sort_counts(self.meta_triples),
        }


def generate_metagraph_streaming(
    nodes_file: Path,
    edges_file: Path,
    graph_name: str,
    graph_version: str | None,
    source_versions: dict[str, str | None] | None = None,
    biolink_version: str | None = None,
) -> MetagraphStats:
    """Generate metagraph statistics from JSONL files using streaming"""
    logging.info(f"Generating metagraph for {graph_name}")

    stats = MetagraphStats(graph_name, graph_version, source_versions, biolink_version)

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

            # Collect edge metadata, keeping each provenance role separate
            if EDGE_PRIMARY_KS in edge:
                stats.primary_knowledge_sources[edge[EDGE_PRIMARY_KS]] += 1
            if EDGE_AGGREGATOR_KS in edge:
                for aggregator_source in edge[EDGE_AGGREGATOR_KS]:
                    stats.aggregator_knowledge_sources[aggregator_source] += 1
            if EDGE_SUPPORTING_SOURCES in edge:
                for supporting_source in edge[EDGE_SUPPORTING_SOURCES]:
                    stats.supporting_data_sources[supporting_source] += 1
            if EDGE_KNOWLEDGE_LEVEL in edge:
                stats.knowledge_levels[edge[EDGE_KNOWLEDGE_LEVEL]] += 1
            if EDGE_AGENT_TYPE in edge:
                stats.agent_types[edge[EDGE_AGENT_TYPE]] += 1
            # Joint knowledge-level x agent-type distribution (for the KLAT heatmap);
            # missing roles fall back to "not_provided" so every edge lands in one cell.
            stats.klat_joint[
                (edge.get(EDGE_KNOWLEDGE_LEVEL, "not_provided"), edge.get(EDGE_AGENT_TYPE, "not_provided"))
            ] += 1

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
    nodes_path: Path,
    edges_path: Path,
    output_dir: Path,
    graph_name: str,
    graph_version: str = None,
    source_versions: dict[str, str | None] | None = None,
    biolink_version: str | None = None,
) -> Path:
    logging.info(f"Generating metagraph for source {graph_name}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate core statistics
    stats = generate_metagraph_streaming(
        nodes_path, edges_path, graph_name, graph_version, source_versions, biolink_version
    )

    # Save main JSON statistics
    file_name = f"{graph_name}_metagraph_{graph_version}.json" if graph_version else f"{graph_name}_metagraph.json"
    json_file_path = output_dir / file_name
    save_metagraph(stats, json_file_path)

    logging.info(f"Metagraph generated for {graph_name}: {json_file_path}")
    return json_file_path
