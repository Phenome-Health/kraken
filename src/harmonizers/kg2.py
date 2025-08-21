"""
RTX-KG2 harmonizer - converts RTX-KG2 format to our schema
"""

from pathlib import Path
import jsonlines
import logging
from ..utils.metagraph import generate_metagraph_for_source

KG2_IGNORE_PROPS = {'domain_range_exclusion'}
KG2_NODE_PROP_NAME_OVERRIDES = {
    "all_categories": "categories",
    "equivalent_curies": "equivalent_ids",
    "all_names": "synonyms",
    "category": "canonical_category"
}
KG2_EDGE_PROP_NAME_OVERRIDES = {
    "id": "kg2c_id",
    "kg2_ids": "kg2pre_ids"
}


def harmonize_kg2(nodes_input: Path, edges_input: Path, nodes_output: Path, edges_output: Path, rules: dict):
    """Harmonize RTX-KG2 to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing RTX-KG2: {nodes_input}, {edges_input} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0

    # Stream and harmonize nodes
    with jsonlines.open(nodes_input, 'r') as reader, jsonlines.open(nodes_output, 'w') as writer:
        for line_num, node in enumerate(reader, 1):
            try:
                harmonized_node = harmonize_kg2_node(node)
                writer.write(harmonized_node)
                node_count += 1
                
                if node_count % 10000 == 0:
                    logging.info(f"Processed {node_count} nodes")
                    
            except (KeyError, TypeError) as e:
                logging.warning(f"Skipping invalid node at line {line_num}: {e}")

    # Stream and harmonize edges
    with jsonlines.open(edges_input, 'r') as reader, jsonlines.open(edges_output, 'w') as writer:
        for line_num, edge in enumerate(reader, 1):
            try:
                # Exclude semmeddb edges and edges with conflicting domain/range
                if not edge.get('domain_range_exclusion') and edge['primary_knowledge_source'] != 'infores:semmeddb':
                    harmonized_edge = harmonize_kg2_edge(edge)
                    writer.write(harmonized_edge)
                    edge_count += 1
                    
                    if edge_count % 10000 == 0:
                        logging.info(f"Processed {edge_count} edges")
                    
            except (KeyError, TypeError) as e:
                logging.warning(f"Skipping invalid edge at line {line_num}: {e}")

    logging.info(f"RTX-KG2 harmonization complete: {node_count} nodes, {edge_count} edges")
    
    # Generate metagraph for harmonized output
    if rules.get('generate_metagraph', True):
        # Store metagraphs in artifacts/metagraphs/harmonized/source_name/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / "kg2"
        
        metagraph_config = rules.get('metagraph_config', {
            'generate_summaries': True,
            'generate_cytoscape': True,
            'generate_html_viewer': True,
            'cytoscape_thresholds': [1, 5, 10]
        })
        
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, "kg2", metagraph_config)
        logging.info("RTX-KG2 metagraph generated")


def harmonize_kg2_node(node: dict) -> dict:
    """Harmonize a single RTX-KG2 node"""
    harmonized_node = {KG2_NODE_PROP_NAME_OVERRIDES.get(property_name, property_name): value
                       for property_name, value in node.items() if property_name not in KG2_IGNORE_PROPS}
    harmonized_node['categories'] = [ensure_biolink_category(category) for category in harmonized_node['categories']]
    harmonized_node['canonical_category'] = ensure_biolink_category(harmonized_node['canonical_category'])
    harmonized_node['provided_by'] = ['infores:kg2']
    return harmonized_node


def harmonize_kg2_edge(edge: dict) -> dict:
    """Harmonize a single RTX-KG2 edge"""
    harmonized_edge = {KG2_EDGE_PROP_NAME_OVERRIDES.get(property_name, property_name): value
                       for property_name, value in edge.items() if property_name not in KG2_IGNORE_PROPS}
    harmonized_edge['aggregator_knowledge_source'] = 'infores:rtx-kg2'
    harmonized_edge['predicate'] = ensure_biolink_predicate(harmonized_edge['predicate'])
    return harmonized_edge


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