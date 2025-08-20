"""
Entity resolution and graph integration functions
"""

from pathlib import Path
from typing import Dict, Optional, List
import logging

from ..utils.kg_io import load_kg, save_kg


def integrate_sources(harmonized_paths: Dict[str, Path], config: dict) -> Path:
    """Merge harmonized sources into unified KG with entity resolution"""
    output_path = Path(config['output_path'])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logging.info("Starting source integration with entity resolution...")

    # Start with RTX-KG2 as the base (highest priority)
    base_kg = load_kg(harmonized_paths['rtx-kg2'])
    logging.info(f"Loaded base KG (RTX-KG2): {base_kg.number_of_nodes()} nodes, {base_kg.number_of_edges()} edges")

    # Build initial node index for entity resolution
    node_index = build_node_index(base_kg)
    logging.info(f"Built initial node index: {len(node_index)} identifiers")

    # Merge each additional source
    for source_name, kg_path in harmonized_paths.items():
        if source_name == 'rtx-kg2':
            continue

        logging.info(f"Merging {source_name}...")
        source_kg = load_kg(kg_path)

        base_kg, node_index = merge_source_with_resolution(
            base_kg, source_kg, node_index, source_name, config
        )

        logging.info(
            f"After merging {source_name}: {base_kg.number_of_nodes()} nodes, {base_kg.number_of_edges()} edges")

    # Save final unified KG
    save_kg(base_kg, output_path)
    logging.info(f"Integration complete! Saved PhenomeKG to: {output_path}")

    return output_path


def build_node_index(kg) -> Dict[str, str]:
    """Build index mapping any identifier to canonical node ID"""
    index = {}

    for node_id, node_data in kg.nodes(data=True):
        equiv_ids = node_data.get('equivalent_ids', [])

        # Map primary ID and all equivalent IDs to this canonical node
        index[node_id] = node_id
        for equiv_id in equiv_ids:
            if equiv_id in index and index[equiv_id] != node_id:
                logging.warning(f"Identifier collision: {equiv_id} maps to both {index[equiv_id]} and {node_id}")
            index[equiv_id] = node_id

    return index


def merge_source_with_resolution(base_kg, source_kg, node_index, source_name, config):
    """Merge source KG into base with entity resolution"""
    namespace_priority = config.get('entity_resolution', {}).get('namespace_priority', [])

    nodes_merged = 0
    nodes_added = 0

    # Process nodes
    for source_node_id, source_node_data in source_kg.nodes(data=True):
        source_equiv_ids = source_node_data.get('equivalent_ids', [])

        # Check if this node matches any existing node
        canonical_id = find_matching_node(source_node_id, source_equiv_ids, node_index)

        if canonical_id:
            # Merge with existing node
            existing_node_data = base_kg.nodes[canonical_id]
            merged_node_data = merge_node_data(
                existing_node_data, source_node_data, source_name, namespace_priority
            )
            base_kg.nodes[canonical_id].update(merged_node_data)

            # Update index with any new equivalent IDs
            for equiv_id in source_equiv_ids:
                if equiv_id not in node_index:
                    node_index[equiv_id] = canonical_id

            nodes_merged += 1
        else:
            # Add as new node
            base_kg.add_node(source_node_id, **source_node_data)

            # Update index
            node_index[source_node_id] = source_node_id
            for equiv_id in source_equiv_ids:
                node_index[equiv_id] = source_node_id

            nodes_added += 1

    # Process edges (with node ID resolution)
    edges_merged = 0
    edges_added = 0

    for source_edge in source_kg.edges(data=True):
        source_subject, source_object, edge_data = source_edge

        # Resolve subject and object to canonical IDs
        canonical_subject = node_index.get(source_subject, source_subject)
        canonical_object = node_index.get(source_object, source_object)

        # Check if edge already exists
        if base_kg.has_edge(canonical_subject, canonical_object):
            # Merge edge data
            existing_edge_data = base_kg.edges[canonical_subject, canonical_object]
            merged_edge_data = merge_edge_data(existing_edge_data, edge_data, source_name)
            base_kg.edges[canonical_subject, canonical_object].update(merged_edge_data)
            edges_merged += 1
        else:
            # Add new edge
            base_kg.add_edge(canonical_subject, canonical_object, **edge_data)
            edges_added += 1

    logging.info(f"  Nodes: {nodes_merged} merged, {nodes_added} added")
    logging.info(f"  Edges: {edges_merged} merged, {edges_added} added")

    return base_kg, node_index


def find_matching_node(node_id: str, equiv_ids: List[str], node_index: Dict[str, str]) -> Optional[str]:
    """Find if this node matches any existing node via equivalent IDs"""
    all_ids = [node_id] + equiv_ids

    for id_val in all_ids:
        if id_val in node_index:
            return node_index[id_val]

    return None


def merge_node_data(existing_node: dict, new_node: dict, source_name: str, namespace_priority: List[str]) -> dict:
    """Merge data from new node into existing node"""
    merged = existing_node.copy()

    # Merge equivalent_ids
    existing_equivs = set(merged.get('equivalent_ids', []))
    new_equivs = set(new_node.get('equivalent_ids', []))
    merged['equivalent_ids'] = list(existing_equivs | new_equivs)

    # Add source provenance
    sources = set(merged.get('sources', [existing_node.get('source', '')]))
    sources.add(source_name)
    merged['sources'] = list(sources)

    # Merge other properties (handle conflicts)
    for key, value in new_node.items():
        if key == 'id':
            # Keep existing canonical ID
            continue
        elif key not in merged or merged[key] is None:
            # Add new property
            merged[key] = value
        elif key not in ['equivalent_ids', 'sources', 'source']:
            # Handle property conflicts
            merged[key] = resolve_property_conflict(
                merged[key], value, key, namespace_priority
            )

    return merged


def merge_edge_data(existing_edge: dict, new_edge: dict, source_name: str) -> dict:
    """Merge edge data, combining sources and handling conflicts"""
    merged = existing_edge.copy()

    # Add source provenance
    sources = set(merged.get('sources', [existing_edge.get('source', '')]))
    sources.add(source_name)
    merged['sources'] = list(sources)

    # For edges, we mostly just want to track multiple sources
    # Could add more sophisticated conflict resolution here if needed

    return merged


def resolve_property_conflict(existing_value, new_value, property_name: str, namespace_priority: List[str]):
    """Resolve conflicts between property values from different sources"""

    # If values are the same, no conflict
    if existing_value == new_value:
        return existing_value

    # For now, just keep the existing value (first source wins)
    # Could implement more sophisticated resolution based on source priority
    logging.debug(f"Property conflict for {property_name}: keeping '{existing_value}' over '{new_value}'")
    return existing_value


def choose_canonical_id(equivalent_ids: set, namespace_priority: List[str] = None) -> str:
    """Choose the canonical ID from a set of equivalent IDs"""
    if not namespace_priority:
        namespace_priority = [
            'CHEMBL.COMPOUND',
            'PUBCHEM.COMPOUND',
            'DRUGBANK',
            'MESH',
            'UNIPROT',
            'ENSEMBL',
            'HGNC'
        ]

    # Try each namespace in priority order
    for namespace in namespace_priority:
        candidates = [id_val for id_val in equivalent_ids if id_val.startswith(namespace)]
        if candidates:
            return sorted(candidates)[0]  # Take first alphabetically if multiple

    # Fallback: shortest ID or lexicographically first
    return min(equivalent_ids, key=lambda x: (len(x), x))
