"""
SPOKE harmonizer - converts SPOKE format to unified Biolink schema
"""

from pathlib import Path
from typing import List
import jsonlines
import logging
from ..utils.metagraph import generate_metagraph_for_source


def harmonize_spoke(input_file: Path, nodes_output: Path, edges_output: Path, rules: dict):
    """Harmonize SPOKE mixed JSONL to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing SPOKE: {input_file} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0

    with jsonlines.open(input_file, 'r') as reader, \
         jsonlines.open(nodes_output, 'w') as nodes_writer, \
         jsonlines.open(edges_output, 'w') as edges_writer:
        
        for line_num, item in enumerate(reader, 1):
            try:
                item_type = item.get('type')
                
                if item_type == 'node':
                    harmonized_node = harmonize_spoke_node(item)
                    nodes_writer.write(harmonized_node)
                    node_count += 1
                    
                    if node_count % 10000 == 0:
                        logging.info(f"Processed {node_count} nodes")
                
                elif item_type == 'relationship':
                    harmonized_edge = harmonize_spoke_edge(item)
                    edges_writer.write(harmonized_edge)
                    edge_count += 1
                    
                    if edge_count % 10000 == 0:
                        logging.info(f"Processed {edge_count} edges")
                        
            except (KeyError, TypeError) as e:
                logging.warning(f"Skipping invalid item at line {line_num}: {e}")

    logging.info(f"SPOKE harmonization complete: {node_count} nodes, {edge_count} edges")
    
    # Generate metagraph for harmonized output
    if rules.get('generate_metagraph', True):
        # Store metagraphs in artifacts/metagraphs/harmonized/source_name/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / "spoke"
        
        metagraph_config = rules.get('metagraph_config', {
            'generate_summaries': True,
            'generate_cytoscape': True,
            'generate_html_viewer': True,
            'cytoscape_thresholds': [1, 5, 10]
        })
        
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, "spoke", metagraph_config)
        logging.info("SPOKE metagraph generated")


def harmonize_spoke_node(node_item: dict) -> dict:
    """Harmonize a single SPOKE node"""
    properties = node_item.get('properties', {})
    labels = node_item.get('labels', [])
    if not labels:
        raise ValueError(f"SPOKE node is missing labels: {node_item}")
    
    # TODO: Stuff other properties into standardized ones..
    harmonized_node = {
        'id': node_item['id'],  # TODO: Convert to standard curies here..
        'categories': map_spoke_labels_to_biolink(labels),
        'name': properties.get('name'),
        'provided_by': ['infores:spoke'],
        'equivalent_ids': [node_item['id']],  # TODO: load any equivalent ids as well..
        'spoke_node': node_item
    }
    return harmonized_node


def harmonize_spoke_edge(edge_item: dict) -> dict:
    """Harmonize a single SPOKE edge"""
    edge_type = edge_item.get('label')
    if not edge_type:
        raise ValueError(f"SPOKE edge is missing type: {edge_item}")
    
     # TODO: Stuff other properties into standardized ones..
    harmonized_edge = {
        'subject': edge_item['start']['id'],  # TODO: Use standard curies here..
        'object': edge_item['start']['id'],
        'predicate': map_spoke_edge_type_to_biolink(edge_type),
        'primary_knowledge_source': "TODO",  # TODO: Replace this placeholder.. 
        'aggregator_knowledge_source': 'infores:spoke',
        'spoke_edge': edge_item
    }
    return harmonized_edge


def map_spoke_labels_to_biolink(labels: List[str]) -> List[str]:
    """Map SPOKE node labels to Biolink categories"""

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
    
    return [label_mapping.get(label, label) for label in labels]


def map_spoke_edge_type_to_biolink(edge_type: str) -> str:
    """Map SPOKE edge types to Biolink predicates"""
    core_edge_type = edge_type.split('_')[0]

    # Simple mapping - extend as needed
    type_mapping = {
        'INTERACTS_WITH': 'biolink:interacts_with',
        'REGULATES': 'biolink:regulates',  # TODO: is this real? think may be old..
        'PART_OF': 'biolink:part_of',
        'TREATS': 'biolink:treats',
        'CAUSES': 'biolink:causes',
        'ASSOCIATED_WITH': 'biolink:associated_with',
        'UPREGULATES': 'biolink:affects'  # TODO: Use qualifiers here... 
    }
    
    return type_mapping.get(core_edge_type, core_edge_type)