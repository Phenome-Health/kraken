"""
SPOKE harmonizer - converts SPOKE format to unified Biolink schema
"""

from pathlib import Path
import json
import logging
from ..utils.metagraph import generate_metagraph_for_source


def harmonize_spoke(input_file: Path, nodes_output: Path, edges_output: Path, rules: dict):
    """Harmonize SPOKE mixed JSONL to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing SPOKE: {input_file} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0

    with open(input_file, 'r') as infile, \
         open(nodes_output, 'w') as nodes_out, \
         open(edges_output, 'w') as edges_out:
        
        for line_num, line in enumerate(infile, 1):
            try:
                item = json.loads(line.strip())
                item_type = item.get('type')
                
                if item_type == 'node':
                    harmonized_node = harmonize_spoke_node(item)
                    nodes_out.write(json.dumps(harmonized_node) + '\n')
                    node_count += 1
                    
                    if node_count % 10000 == 0:
                        logging.info(f"Processed {node_count} nodes")
                
                elif item_type == 'edge':
                    harmonized_edge = harmonize_spoke_edge(item)
                    edges_out.write(json.dumps(harmonized_edge) + '\n')
                    edge_count += 1
                    
                    if edge_count % 10000 == 0:
                        logging.info(f"Processed {edge_count} edges")
                        
            except (json.JSONDecodeError, KeyError) as e:
                logging.warning(f"Skipping invalid item at line {line_num}: {e}")

    logging.info(f"SPOKE harmonization complete: {node_count} nodes, {edge_count} edges")
    
    # Generate metagraph for harmonized output
    if rules.get('generate_metagraph', True):
        metagraph_dir = nodes_output.parent / "metagraphs"
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, "spoke")
        logging.info("SPOKE metagraph generated")


def harmonize_spoke_node(node_item: dict) -> dict:
    """Harmonize a single SPOKE node"""
    properties = node_item.get('properties', {})
    labels = node_item.get('labels', [])
    
    # Map SPOKE labels to Biolink categories
    biolink_category = map_spoke_labels_to_biolink(labels)
    
    harmonized = {
        'id': properties.get('identifier') or node_item.get('id'),
        'category': biolink_category,
        'name': properties.get('name'),
        'source': 'spoke'
    }
    
    # Copy other properties
    for key, value in properties.items():
        if key not in ['identifier', 'name'] and key not in harmonized:
            harmonized[key] = value
    
    # Add original labels for reference
    if labels:
        harmonized['spoke_labels'] = labels
    
    return harmonized


def harmonize_spoke_edge(edge_item: dict) -> dict:
    """Harmonize a single SPOKE edge"""
    properties = edge_item.get('properties', {})
    
    harmonized = {
        'subject': edge_item.get('startNode') or edge_item.get('source'),
        'object': edge_item.get('endNode') or edge_item.get('target'),  
        'predicate': map_spoke_edge_type_to_biolink(edge_item.get('type')),
        'source': 'spoke'
    }
    
    # Copy other properties
    for key, value in properties.items():
        if key not in harmonized:
            harmonized[key] = value
    
    return harmonized


def map_spoke_labels_to_biolink(labels: list) -> str:
    """Map SPOKE node labels to Biolink categories"""
    if not labels:
        return "biolink:NamedThing"
    
    # Simple mapping - extend as needed
    label_mapping = {
        'Anatomy': 'biolink:AnatomicalEntity',
        'Gene': 'biolink:Gene', 
        'Protein': 'biolink:Protein',
        'Disease': 'biolink:Disease',
        'Compound': 'biolink:ChemicalEntity',
        'Pathway': 'biolink:Pathway',
        'BiologicalProcess': 'biolink:BiologicalProcess',
        'MolecularFunction': 'biolink:MolecularActivity',
        'CellularComponent': 'biolink:CellularComponent'
    }
    
    # Use first recognized label
    for label in labels:
        if label in label_mapping:
            return label_mapping[label]
    
    return "biolink:NamedThing"


def map_spoke_edge_type_to_biolink(edge_type: str) -> str:
    """Map SPOKE edge types to Biolink predicates"""
    if not edge_type:
        return "biolink:related_to"
    
    # Simple mapping - extend as needed
    type_mapping = {
        'INTERACTS_WITH': 'biolink:interacts_with',
        'REGULATES': 'biolink:regulates', 
        'PART_OF': 'biolink:part_of',
        'TREATS': 'biolink:treats',
        'CAUSES': 'biolink:causes',
        'ASSOCIATED_WITH': 'biolink:associated_with'
    }
    
    return type_mapping.get(edge_type, f"biolink:{edge_type.lower()}")