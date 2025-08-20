"""
Biomapper export utilities
Exports node type-specific files for biomapper module
"""

from pathlib import Path
import logging
from collections import defaultdict
from ..utils.kg_io import save_kg, create_kg_from_nodes_edges


def export_for_biomapper(unified_kg: nx.MultiDiGraph, config: dict):
    """Export unified KG divided by node types for biomapper"""
    output_dir = Path(config['output_directory'])
    output_dir.mkdir(parents=True, exist_ok=True)

    node_types = config.get('node_types_to_export', [])
    include_edges = config.get('include_edges', True)

    logging.info(f"Exporting {len(node_types)} node types for biomapper...")

    # Group nodes by type
    nodes_by_type = group_nodes_by_type(unified_kg, node_types)

    # Export each node type as a separate file
    for node_type, nodes in nodes_by_type.items():
        if nodes:  # Only export if we have nodes of this type
            export_node_type_subgraph(
                unified_kg, nodes, node_type, output_dir, include_edges
            )

    # Also create a summary file
    create_biomapper_summary(nodes_by_type, output_dir)

    logging.info(f"Biomapper export complete: {len(nodes_by_type)} files in {output_dir}")


def group_nodes_by_type(kg: nx.MultiDiGraph, target_types: list) -> dict:
    """Group nodes by their biolink category"""
    nodes_by_type = defaultdict(list)

    for node_id, node_data in kg.nodes(data=True):
        category = node_data.get('category', 'biolink:NamedThing')

        # Handle both single categories and lists
        if isinstance(category, list):
            categories = category
        else:
            categories = [category]

        # Check if any of the node's categories are in our target list
        for cat in categories:
            if cat in target_types:
                nodes_by_type[cat].append(node_id)
                break  # Only add to one category to avoid duplicates

    return dict(nodes_by_type)


def export_node_type_subgraph(kg: nx.MultiDiGraph, node_ids: list, node_type: str,
                              output_dir: Path, include_edges: bool):
    """Export a subgraph containing only nodes of a specific type"""

    # Create subgraph with these nodes
    subgraph_nodes = []
    for node_id in node_ids:
        node_data = kg.nodes[node_id]
        subgraph_nodes.append({'id': node_id, **node_data})

    subgraph_edges = []
    if include_edges:
        # Include edges between nodes of this type, or edges connecting to other types
        for node_id in node_ids:
            # Get all edges involving this node
            for neighbor in kg.neighbors(node_id):
                edge_data = kg.edges[node_id, neighbor]
                subgraph_edges.append({
                    'subject': node_id,
                    'object': neighbor,
                    **edge_data
                })

            # Also get incoming edges
            for predecessor in kg.predecessors(node_id):
                if kg.has_edge(predecessor, node_id):
                    edge_data = kg.edges[predecessor, node_id]
                    subgraph_edges.append({
                        'subject': predecessor,
                        'object': node_id,
                        **edge_data
                    })

    # Create the subgraph
    subgraph = create_kg_from_nodes_edges(subgraph_nodes, subgraph_edges)

    # Save with a clean filename
    clean_type_name = node_type.replace('biolink:', '').replace(':', '_')
    output_path = output_dir / f"{clean_type_name}.json"

    save_kg(subgraph, output_path)
    logging.info(f"  Exported {len(subgraph_nodes)} {node_type} nodes to {output_path}")


def create_biomapper_summary(nodes_by_type: dict, output_dir: Path):
    """Create a summary file with statistics about the export"""
    summary = {
        'export_summary': {
            'total_node_types': len(nodes_by_type),
            'node_counts_by_type': {
                node_type: len(nodes)
                for node_type, nodes in nodes_by_type.items()
            },
            'total_nodes_exported': sum(len(nodes) for nodes in nodes_by_type.values())
        }
    }

    import json
    summary_path = output_dir / "export_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logging.info(f"Export summary saved to: {summary_path}")