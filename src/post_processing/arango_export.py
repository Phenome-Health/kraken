"""
ArangoDB export utilities
Prepares the unified KG for import into ArangoDB
"""

from pathlib import Path
import re
import logging
from ..utils.kg_io import save_kg


def prepare_for_arango(unified_kg: nx.MultiDiGraph, config: dict):
    """Prepare unified KG for ArangoDB import"""
    output_path = Path(config['output_path'])
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logging.info("Preparing KG for ArangoDB...")

    # Create a copy to modify
    arango_kg = unified_kg.copy()

    # Add _key fields to all nodes
    add_arango_keys_to_nodes(arango_kg, config.get('key_field_config', {}))

    # Add _key fields to all edges (if needed)
    add_arango_keys_to_edges(arango_kg, config.get('key_field_config', {}))

    # Save the ArangoDB-ready version
    save_kg(arango_kg, output_path)

    logging.info(f"ArangoDB export saved to: {output_path}")


def add_arango_keys_to_nodes(kg: nx.MultiDiGraph, key_config: dict):
    """Add _key field to all nodes according to ArangoDB requirements"""
    remove_invalid = key_config.get('remove_invalid_chars', True)
    max_length = key_config.get('max_length', 254)
    prefix_numbers = key_config.get('prefix_numbers', True)

    for node_id, node_data in kg.nodes(data=True):
        arango_key = create_arango_key(node_id, remove_invalid, max_length, prefix_numbers)
        node_data['_key'] = arango_key


def add_arango_keys_to_edges(kg: nx.MultiDiGraph, key_config: dict):
    """Add _key field to all edges if needed"""
    # ArangoDB can auto-generate edge keys, but you might want custom ones
    for u, v, edge_data in kg.edges(data=True):
        # Create a deterministic key based on subject, object, and predicate
        predicate = edge_data.get('predicate', 'related_to')
        edge_key = f"{u}_{predicate}_{v}"
        arango_key = create_arango_key(edge_key, True, 254, True)
        edge_data['_key'] = arango_key


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