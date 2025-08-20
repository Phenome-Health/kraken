"""
RTX-KG2 harmonizer - converts RTX-KG2 format to our schema
"""

from pathlib import Path
import logging
from ..utils.kg_io import load_kg, save_kg, stream_kg_nodes, create_kg_from_nodes_edges


def harmonize_kg2(input_path: Path, output_path: Path, rules: dict):
    """Harmonize RTX-KG2 to unified Biolink schema"""
    logging.info(f"Harmonizing RTX-KG2: {input_path} -> {output_path}")

    # RTX-KG2 might already be in good Biolink format, so this might be minimal
    # Just standardize any format differences

    harmonized_nodes = []
    harmonized_edges = []

    kg = load_kg(input_path)

    for node in kg.nodes(data=True):
        node_id, node_data = node
        harmonized_node = harmonize_kg2_node(node_id, node_data)
        harmonized_nodes.append(harmonized_node)

    for edge in kg.edges(data=True):
        source, target, edge_data = edge
        harmonized_edge = harmonize_kg2_edge(source, target, edge_data)
        harmonized_edges.append(harmonized_edge)

    # Save harmonized graph
    harmonized_kg = create_kg_from_nodes_edges(harmonized_nodes, harmonized_edges)
    save_kg(harmonized_kg, output_path)

    logging.info(f"RTX-KG2 harmonization complete: {len(harmonized_nodes)} nodes, {len(harmonized_edges)} edges")


def harmonize_kg2_node(node_id: str, node_data: dict) -> dict:
    """Harmonize a single RTX-KG2 node"""
    harmonized = {
        'id': node_id,
        'category': ensure_biolink_category(node_data.get('category')),
        'name': node_data.get('name'),
        'equivalent_identifiers': node_data.get('equivalent_identifiers', []),
        'source': 'rtx-kg2'
    }

    # Copy other properties
    for key, value in node_data.items():
        if key not in harmonized:
            harmonized[key] = value

    return harmonized


def harmonize_kg2_edge(source: str, target: str, edge_data: dict) -> dict:
    """Harmonize a single RTX-KG2 edge"""
    return {
        'subject': source,
        'object': target,
        'predicate': ensure_biolink_predicate(edge_data.get('predicate')),
        'source': 'rtx-kg2',
        **{k: v for k, v in edge_data.items() if k not in ['subject', 'object', 'predicate']}
    }


def ensure_biolink_category(category):
    """Ensure category is in proper biolink format"""
    if not category:
        return "biolink:NamedThing"

    if not category.startswith("biolink:"):
        return f"biolink:{category}"

    return category


def ensure_biolink_predicate(predicate):
    """Ensure predicate is in proper biolink format"""
    if not predicate:
        return "biolink:related_to"

    if not predicate.startswith("biolink:"):
        return f"biolink:{predicate}"

    return predicate