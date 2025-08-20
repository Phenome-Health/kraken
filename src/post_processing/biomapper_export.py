"""
Biomapper export utilities
Exports node type-specific files for biomapper module
"""

from pathlib import Path
import logging
import json
from collections import defaultdict
from ..utils.kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl


def export_for_biomapper_streaming(nodes_file: Path, edges_file: Path, output_dir: Path, config: dict):
    """Export unified KG divided by node types for biomapper using streaming"""
    output_dir.mkdir(parents=True, exist_ok=True)

    node_types = config.get('node_types_to_export', [])
    include_edges = config.get('include_edges', True)

    logging.info(f"Streaming export of {len(node_types)} node types for biomapper...")

    # First pass: group nodes by type and collect node IDs for each type
    nodes_by_type = defaultdict(list)
    node_type_membership = {}  # node_id -> node_type (for edge processing)
    
    for node in stream_nodes_from_jsonl(nodes_file):
        node_id = node['id']
        category = node.get('category', 'biolink:NamedThing')

        # Handle both single categories and lists
        if isinstance(category, list):
            categories = category
        else:
            categories = [category]

        # Check if any of the node's categories are in our target list
        for cat in categories:
            if cat in node_types:
                nodes_by_type[cat].append(node)
                node_type_membership[node_id] = cat
                break  # Only add to one category to avoid duplicates

    logging.info(f"Found nodes in {len(nodes_by_type)} categories")

    # Export each node type as a separate file
    for node_type, nodes in nodes_by_type.items():
        if nodes:  # Only export if we have nodes of this type
            export_node_type_files(nodes, node_type, output_dir, edges_file, 
                                  node_type_membership, include_edges)

    # Create a summary file
    create_biomapper_summary({k: len(v) for k, v in nodes_by_type.items()}, output_dir)

    logging.info(f"Biomapper export complete: {len(nodes_by_type)} files in {output_dir}")


def export_node_type_files(nodes: list, node_type: str, output_dir: Path, 
                          edges_file: Path, node_type_membership: dict, include_edges: bool):
    """Export nodes of a specific type and their associated edges"""
    
    # Clean filename
    clean_type_name = node_type.replace('biolink:', '').replace(':', '_')
    nodes_output = output_dir / f"{clean_type_name}_nodes.jsonl"
    edges_output = output_dir / f"{clean_type_name}_edges.jsonl"
    
    # Save nodes
    with open(nodes_output, 'w') as f:
        for node in nodes:
            f.write(json.dumps(node) + '\n')
    
    logging.info(f"  Exported {len(nodes)} {node_type} nodes to {nodes_output}")
    
    if include_edges:
        # Collect node IDs for this type
        node_ids_in_type = {node['id'] for node in nodes}
        
        # Stream through edges and save those involving nodes of this type
        edge_count = 0
        with open(edges_output, 'w') as f:
            for edge in stream_edges_from_jsonl(edges_file):
                subject = edge.get('subject')
                obj = edge.get('object')
                
                # Include edge if either endpoint is in this node type
                if subject in node_ids_in_type or obj in node_ids_in_type:
                    f.write(json.dumps(edge) + '\n')
                    edge_count += 1
        
        logging.info(f"  Exported {edge_count} edges to {edges_output}")


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


def create_biomapper_summary(node_counts_by_type: dict, output_dir: Path):
    """Create a summary file with statistics about the export"""
    summary = {
        'export_summary': {
            'total_node_types': len(node_counts_by_type),
            'node_counts_by_type': node_counts_by_type,
            'total_nodes_exported': sum(node_counts_by_type.values())
        }
    }

    summary_path = output_dir / "export_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logging.info(f"Export summary saved to: {summary_path}")