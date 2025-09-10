"""
RTX-KG2 harmonizer - converts RTX-KG2 format to our schema
"""

from pathlib import Path
import jsonlines
import logging
from ..utils.metagraph import generate_metagraph_for_source
from ..utils.constants import *
from ..utils.kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl

KG2_IGNORE_PROPS = {'domain_range_exclusion'}
KG2_NODE_PROP_NAME_OVERRIDES = {
    "all_categories": CATEGORIES,
    "equivalent_curies": EQUIVALENT_IDS,
    "all_names": SYNONYMS,
    "category": "canonical_category"
}
KG2_EDGE_PROP_NAME_OVERRIDES = {
    "id": "kg2c_id",
    "kg2_ids": "kg2pre_ids"
}


def harmonize_kg2(nodes_input: Path, edges_input: Path, nodes_output: Path, edges_output: Path, biolink_version: str, build_metagraph: bool):
    """Harmonize RTX-KG2 to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing RTX-KG2: {nodes_input}, {edges_input} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0

    # Stream and harmonize nodes
    with jsonlines.open(nodes_output, 'w') as writer:
        for node in stream_nodes_from_jsonl(nodes_input):
            harmonized_node = harmonize_kg2_node(node)
            writer.write(harmonized_node)
            node_count += 1

    # Stream and harmonize edges
    with jsonlines.open(edges_output, 'w') as writer:
        for edge in stream_edges_from_jsonl(edges_input):
            # Exclude semmeddb edges and edges with conflicting domain/range
            if not edge.get('domain_range_exclusion') and edge['primary_knowledge_source'] != 'infores:semmeddb':
                harmonized_edge = harmonize_kg2_edge(edge)
                writer.write(harmonized_edge)
                edge_count += 1

    logging.info(f"RTX-KG2 harmonization complete: {node_count} nodes, {edge_count} edges")
    
    if build_metagraph:
        # Generate metagraph for harmonized output
        # Store metagraphs in artifacts/metagraphs/harmonized/source_name/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / "kg2"
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, "kg2")
        logging.info("RTX-KG2 metagraph generated")


def harmonize_kg2_node(node: dict) -> dict:
    """Harmonize a single RTX-KG2 node"""
    harmonized_node = {KG2_NODE_PROP_NAME_OVERRIDES.get(property_name, property_name): value
                       for property_name, value in node.items() if property_name not in KG2_IGNORE_PROPS}
    harmonized_node[PROVIDED_BY] = [KG2_INFORES]

    # Handle missing equivalent_curies property (happens with KG2pre build node)
    if EQUIVALENT_IDS not in harmonized_node:
        harmonized_node[EQUIVALENT_IDS] = [harmonized_node[ID]]
    harmonized_node[EQUIVALENT_IDS] = sorted(harmonized_node[EQUIVALENT_IDS], key=str.casefold)

    # Ensure all nodes have 'synonyms' as expected, and 'name' is in it
    if harmonized_node.get(NAME):
        if not harmonized_node.get(SYNONYMS):
            harmonized_node[SYNONYMS] = [harmonized_node[NAME]]
        if harmonized_node[NAME] not in harmonized_node[SYNONYMS]:
            harmonized_node[SYNONYMS].append(harmonized_node[NAME])

    return harmonized_node


def harmonize_kg2_edge(edge: dict) -> dict:
    """Harmonize a single RTX-KG2 edge"""
    harmonized_edge = {KG2_EDGE_PROP_NAME_OVERRIDES.get(property_name, property_name): value
                       for property_name, value in edge.items() if property_name not in KG2_IGNORE_PROPS}
    harmonized_edge[AGGREGATOR_KS] = KG2_INFORES

    # Make the KG2c edge ID property into a list, because we want it merged during entity resolution (if necessary)
    harmonized_edge['kg2c_ids'] = [harmonized_edge['kg2c_id']]
    del harmonized_edge['kg2c_id']

    return harmonized_edge
