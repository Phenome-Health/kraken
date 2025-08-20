"""
Knowledge graph I/O utilities
"""

from pathlib import Path
import json
import logging


def load_kg(kg_path: Path) -> nx.MultiDiGraph:
    """Load knowledge graph from file"""
    logging.debug(f"Loading KG from {kg_path}")

    if kg_path.suffix == '.json':
        return load_kg_from_json(kg_path)
    else:
        raise ValueError(f"Unsupported KG format: {kg_path.suffix}")


def save_kg(kg: nx.MultiDiGraph, output_path: Path):
    """Save knowledge graph to file"""
    logging.debug(f"Saving KG to {output_path}")

    if output_path.suffix == '.json':
        save_kg_to_json(kg, output_path)
    else:
        raise ValueError(f"Unsupported output format: {output_path.suffix}")


def load_kg_from_json(json_path: Path) -> nx.MultiDiGraph:
    """Load KG from JSON format"""
    with open(json_path) as f:
        data = json.load(f)

    kg = nx.MultiDiGraph()

    # Add nodes
    if 'nodes' in data:
        for node in data['nodes']:
            node_id = node['id']
            kg.add_node(node_id, **{k: v for k, v in node.items() if k != 'id'})

    # Add edges
    if 'edges' in data:
        for edge in data['edges']:
            subject = edge['subject']
            obj = edge['object']
            kg.add_edge(subject, obj, **{k: v for k, v in edge.items() if k not in ['subject', 'object']})

    return kg


def save_kg_to_json(kg: nx.MultiDiGraph, json_path: Path):
    """Save KG to JSON format"""

    # Convert to JSON-serializable format
    nodes = []
    for node_id, node_data in kg.nodes(data=True):
        nodes.append({'id': node_id, **node_data})

    edges = []
    for subject, obj, edge_data in kg.edges(data=True):
        edges.append({'subject': subject, 'object': obj, **edge_data})

    data = {
        'nodes': nodes,
        'edges': edges
    }

    with open(json_path, 'w') as f:
        json.dump(data, f, indent=2)


def stream_kg_nodes(kg_path: Path):
    """Stream nodes from KG file without loading entire graph into memory"""
    # This would be more complex for very large files
    # For now, just load and iterate
    kg = load_kg(kg_path)
    for node_id, node_data in kg.nodes(data=True):
        yield {'id': node_id, **node_data}


def create_kg_from_nodes_edges(nodes: list, edges: list) -> nx.MultiDiGraph:
    """Create NetworkX graph from lists of nodes and edges"""
    kg = nx.MultiDiGraph()

    for node in nodes:
        node_id = node['id']
        kg.add_node(node_id, **{k: v for k, v in node.items() if k != 'id'})

    for edge in edges:
        kg.add_edge(edge['subject'], edge['object'],
                    **{k: v for k, v in edge.items() if k not in ['subject', 'object']})

    return kg
