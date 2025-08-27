"""
ArangoDB export utilities
Prepares the unified KG for import into ArangoDB
"""

from collections import defaultdict
import html
from pathlib import Path
import re
import logging
from typing import Tuple
import unicodedata
import jsonlines

from ..utils.constants import *
from ..utils.kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl
from ..utils.general import create_edge_key, clean_key_for_arango
from bmt import Toolkit


IGNORE_PROPS = {}
NODE_PROP_NAME_OVERRIDES = {
    "categories": "entity_types",
    "canonical_category": "canonical_entity_type"
}
EDGE_PROP_NAME_OVERRIDES = {
    "predicate": "connection_type"
}


def prepare_for_arango(nodes_path: Path, edges_path: Path, output_dir: Path, config: dict, biolink_version: str) -> Tuple[Path, Path]:
    """Prepare unified KG for ArangoDB import using streaming"""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    arango_nodes_path = output_dir / "arango_nodes.jsonl"
    arango_edges_path = output_dir / "arango_edges.jsonl"
    
    logging.info("Preparing KG for ArangoDB...")
    
    logging.info(f"Initiating Biolink Model Toolkit..")
    biolink_url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{biolink_version}/biolink-model.yaml"
    bmt = Toolkit(schema=biolink_url)
    ancestor_map = defaultdict(set)

    # First tally up neighbor counts
    neighbor_counts, neighbor_counts_by_type = count_neighbors(nodes_path, edges_path)
    
    # Process nodes
    node_count = 0
    with jsonlines.open(arango_nodes_path, 'w') as writer:
        for node in stream_nodes_from_jsonl(nodes_path):
            arango_node = create_arango_node(node, ancestor_map, bmt, neighbor_counts, neighbor_counts_by_type)
            writer.write(arango_node)

    # Process edges
    edge_count = 0
    with jsonlines.open(arango_edges_path, 'w') as writer:
        for edge in stream_edges_from_jsonl(edges_path):
            arango_edge = create_arango_edge(edge, ancestor_map, bmt)
            writer.write(arango_edge)

    logging.info(f"ArangoDB export complete: {node_count} nodes, {edge_count} edges")
    logging.info(f"Files saved to: {arango_nodes_path}, {arango_edges_path}")
    
    return arango_nodes_path, arango_edges_path


def clean_text(text: any) -> str:
    if not isinstance(text, str):
        # Handle weird case where some KG2 nodes have a name of True (a bool)
        text = str(text)
    unescaped_text = html.unescape(text)
    cleaned_text = unicodedata.normalize('NFC', unescaped_text)
    return cleaned_text


def normalize_text(text: str) -> str:
    lowercase_text = text.lower()

    # Replace all non-alphanumeric characters with a space
    # [^a-zA-Z0-9] matches any character that is NOT a letter or digit
    no_punct = re.sub(r'[^a-z0-9\s]', ' ', lowercase_text)

    # Replace multiple spaces with a single space and strip leading/trailing spaces
    cleaned_text = re.sub(r'\s+', ' ', no_punct).strip()

    return cleaned_text


def custom_sort_key(item) -> Tuple[int, str]:
    """
    Sort key that groups items starting with letters first, then items starting with numbers.
    Sorting within each group is case-insensitive.
    """
    s_item = str(item).lower() # Convert to lowercase string for comparison

    if s_item and s_item[0].isalpha():
        # Group 0: Items starting with a letter
        return 0, s_item
    else:
        # Group 1: Items starting with a number/other symbol
        return 1, s_item


def create_arango_node(node: dict, ancestor_map: defaultdict[set], bmt: Toolkit, neighbor_counts: dict, neighbor_counts_by_type: dict) -> dict:
    # Create arango version of node
    arango_node = {NODE_PROP_NAME_OVERRIDES.get(property_name, property_name): value
                    for property_name, value in node.items() if property_name not in IGNORE_PROPS}
    arango_node['_key'] = clean_key_for_arango(node['id'])

    # Clean name and make normalized versions of it
    if arango_node.get('name'):
        arango_node['name'] = clean_text(arango_node['name'])
        arango_node['name_normalized'] = normalize_text(arango_node['name'])
    
    # Clean synonyms and make normalized version of it
    if arango_node.get('synonyms'):
        arango_node['synonyms'] = sorted(list({clean_text(synonym)
                                                for synonym in arango_node['synonyms']}),
                                            key=custom_sort_key)
        arango_node['synonyms_normalized'] = sorted(list({normalize_text(synonym)
                                                            for synonym in arango_node['synonyms']}),
                                                    key=custom_sort_key)

    # Clean the description field, if present (get rid of html-encoded characters)
    if arango_node.get('description'):
        arango_node['description'] = clean_text(arango_node['description'])

    # Add expanded entity types (includes Biolink ancestors)
    ancestors = set()
    for entity_type in arango_node['entity_types']:
        if entity_type not in ancestor_map:
            ancestor_map[entity_type] = set(bmt.get_ancestors(entity_type,
                                                                formatted=True,
                                                                mixin=True,
                                                                reflexive=True))
        ancestors |= ancestor_map[entity_type]
    arango_node["entity_types_ancestral"] = list(ancestors)

    # Extract prefixes from equivalent IDs for easy lookup later
    if 'equivalent_ids' not in arango_node:
        print(arango_node)
    prefixes = {equivalent_id.split(":")[0] for equivalent_id in arango_node["equivalent_ids"]}
    arango_node["id_prefixes"] = sorted(list(prefixes))

    # Annotate nodes with their neighbor counts (calculated during edges conversion)
    arango_node["num_neighbors"] = neighbor_counts.get(arango_node["id"], 0)
    arango_node["neighbor_counts"] = neighbor_counts_by_type.get(arango_node["id"], dict())
    
    return arango_node


def create_arango_edge(edge: dict, ancestor_map: dict, bmt: Toolkit) -> dict:
    arango_edge = {EDGE_PROP_NAME_OVERRIDES.get(prop_name, prop_name): value
                    for prop_name, value in edge.items()}
    arango_edge["_key"] = create_edge_key(edge)
    arango_edge["_from"] = f"nodes/{clean_key_for_arango(edge[SUBJECT])}"
    arango_edge["_to"] = f"nodes/{clean_key_for_arango(edge[OBJECT])}"

    # Add expanded connection types (includes Biolink ancestors)
    connection_type = arango_edge["connection_type"]
    if connection_type not in ancestor_map:
        ancestor_map[connection_type] = set(bmt.get_ancestors(connection_type,
                                                                formatted=True,
                                                                mixin=True,
                                                                reflexive=True))
    arango_edge["connection_types_ancestral"] = list(ancestor_map[connection_type])
    return arango_edge


def count_neighbors(nodes_file_path: str, edges_file_path: str) -> Tuple[defaultdict, defaultdict]:
    neighbor_counts = defaultdict(int)
    neighbor_counts_by_type = defaultdict(lambda: defaultdict(int))

    logging.info(f"Loading categories map")
    with jsonlines.open(nodes_file_path, "r") as reader:
        categories_map = {node[ID]: node[CATEGORIES] for node in reader}

    # Tally up neighbor counts
    logging.info(f"Beginning to count neighbors")
    with jsonlines.open(edges_file_path, "r") as reader:
        for edge in reader:
            subject_id = edge[SUBJECT]
            object_id = edge[OBJECT]
            neighbor_counts[subject_id] += 1
            neighbor_counts[object_id] += 1
            if subject_id in categories_map and object_id in categories_map:
                for subj_category in categories_map[subject_id]:
                    neighbor_counts_by_type[object_id][subj_category] += 1
                for obj_category in categories_map[object_id]:
                    neighbor_counts_by_type[subject_id][obj_category] += 1
            else:
                logging.warning(f"Edge from {subject_id} to {object_id} is an orphan")

    return neighbor_counts, neighbor_counts_by_type
