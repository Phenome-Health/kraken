"""
RTX-KG2 harmonizer - converts RTX-KG2 format to our schema
"""

from pathlib import Path
import jsonlines
import logging
from ..utils.constants import *
from ..utils.kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl
from ..utils.general import create_node, create_edge


def harmonize_kg2(nodes_input: Path, edges_input: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    """Harmonize RTX-KG2 to unified Biolink schema using streaming"""
    logging.info(f"Harmonizing RTX-KG2: {nodes_input}, {edges_input} -> {nodes_output}, {edges_output}")

    node_count = 0
    edge_count = 0

    # Stream and harmonize nodes
    with jsonlines.open(nodes_output, 'w') as writer:
        for node in stream_nodes_from_jsonl(nodes_input):
            harmonized_node = create_node(
                curie=node['id'],
                categories=node['all_categories'],
                provided_by=[KG2_INFORES],
                equivalent_ids=node.get('equivalent_curies', [node['id']]),
                name=node.get('name'),
                synonyms=node.get('all_names'),
                iri=node.get('iri'),
                description=node.get('description'),
                attributes={KG2_INFORES: {
                    'canonical_category': node['category']
                }}
            )
            writer.write(harmonized_node)
            node_count += 1

    # Stream and harmonize edges
    with jsonlines.open(edges_output, 'w') as writer:
        for edge in stream_edges_from_jsonl(edges_input):
            # Exclude semmeddb edges and edges with conflicting domain/range
            if not edge.get('domain_range_exclusion') and edge['primary_knowledge_source'] != 'infores:semmeddb':
                harmonized_edge = create_edge(
                    subject_id=edge['subject'],
                    object_id=edge['object'],
                    predicate=edge['predicate'],
                    primary_ks=edge['primary_knowledge_source'],
                    knowledge_level=edge.get('knowledge_level', 'not_provided'),
                    agent_type=edge.get('agent_type', 'not_provided'),
                    aggregator_ks=KG2_INFORES,
                    qualified_predicate=edge.get('qualified_predicate'),
                    qualified_direction=edge.get('qualified_object_direction'),
                    qualified_aspect=edge.get('qualified_object_aspect'),
                    publications=edge.get('publications'),
                    publications_info = edge.get('publications_info'),
                    attributes={KG2_INFORES: {
                        'kg2c_ids': [edge['id']],
                        'kg2pre_ids': edge['kg2_ids']
                    }}
                )
                writer.write(harmonized_edge)
                edge_count += 1

    logging.info(f"RTX-KG2 harmonization complete: {node_count} nodes, {edge_count} edges")
