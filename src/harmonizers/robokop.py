"""
ROBOKOP harmonizer - converts to our unified schema
"""

from pathlib import Path
import jsonlines
import logging
from ..utils.constants import *
from ..utils.kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl
from ..utils.general import create_node, create_edge

from bmt import Toolkit


def harmonize_robokop(nodes_input: Path, edges_input: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    """Harmonize ROBOKOP to our unified schema using streaming"""
    logging.info(f"Harmonizing ROBOKOP: {nodes_input}, {edges_input} -> {nodes_output}, {edges_output}")

    biolink_url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{biolink_version}/biolink-model.yaml"
    logging.info(f"Initializing bmt (Biolink Model Toolkit) for version {biolink_version}...")
    bmt = Toolkit(schema=biolink_url)

    node_count = 0
    edge_count = 0

    core_robokop_node_props = {'id', 'name', 'category', 'description', 'equivalent_identifiers'}

    # Stream and harmonize nodes
    with jsonlines.open(nodes_output, 'w') as writer:
        for node in stream_nodes_from_jsonl(nodes_input):

            # TODO: check if they have any synonym properties?

            # Remove the purely ancestral categories (we do that expansion at query time)
            all_proper_ancestors = set()
            categories = set(node['category'])
            for category in categories:
                proper_ancestors = set(bmt.get_ancestors(category,
                                                         formatted=True,
                                                         mixin=True,
                                                         reflexive=False))
                all_proper_ancestors |= proper_ancestors
            leaf_categories = categories.difference(all_proper_ancestors)

            harmonized_node = create_node(
                curie=node['id'],
                categories=list(leaf_categories),
                provided_by=[ROBOKOP_INFORES],
                equivalent_ids=node['equivalent_identifiers'] if node.get('equivalent_identifiers') else [node['id']],
                name=node.get('name'),
                description=node.get('description'),
                attributes={ROBOKOP_INFORES: {
                    prop_name: value for prop_name, value in node.items() if prop_name not in core_robokop_node_props
                }}
            )
            writer.write(harmonized_node)
            node_count += 1

    core_robokop_edge_props = CORE_EDGE_PROPERTIES.union({'object_direction_qualifier',
                                                          'object_aspect_qualifier',
                                                          'publications', 'sentences'})

    # Stream and harmonize edges
    with jsonlines.open(edges_output, 'w') as writer:
        for edge in stream_edges_from_jsonl(edges_input):
            harmonized_edge = create_edge(
                subject_id=edge['subject'],
                object_id=edge['object'],
                predicate=edge['predicate'],
                primary_ks=edge['primary_knowledge_source'],
                knowledge_level=edge.get('knowledge_level', 'not_provided'),
                agent_type=edge.get('agent_type', 'not_provided'),
                aggregator_ks=ROBOKOP_INFORES,
                qualified_predicate=edge.get('qualified_predicate'),
                qualified_direction=edge.get('object_direction_qualifier'),
                qualified_aspect=edge.get('object_aspect_qualifier'),
                publications=edge.get('publications'),
                publications_info = edge.get('sentences'),
                attributes={ROBOKOP_INFORES: {
                    prop_name: value for prop_name, value in edge.items() if prop_name not in core_robokop_edge_props
                }}
            )
            writer.write(harmonized_edge)
            edge_count += 1

    logging.info(f"ROBOKOP harmonization complete: {node_count} nodes, {edge_count} edges")
