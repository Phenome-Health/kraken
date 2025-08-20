"""
Metagraph generation utilities for Biolink knowledge graphs
Analyzes node categories, edge predicates, and connectivity patterns
"""

from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, Set, Tuple, Any, Iterator
import json
import logging

from .kg_io import stream_nodes_from_jsonl, stream_edges_from_jsonl


class MetagraphStats:
    """Container for metagraph statistics"""
    
    def __init__(self, source_name: str = "unknown"):
        self.source_name = source_name
        
        # Node statistics
        self.node_categories = Counter()  # category -> count
        self.total_nodes = 0
        
        # Edge statistics  
        self.edge_predicates = Counter()  # predicate -> count
        self.category_pairs = Counter()  # (subject_category, object_category) -> count
        self.predicate_category_pairs = Counter()  # (predicate, subject_cat, object_cat) -> count
        self.total_edges = 0
        
        # Connectivity statistics
        self.node_degrees = defaultdict(int)  # node_id -> degree
        self.category_connectivity = defaultdict(lambda: defaultdict(int))  # subj_cat -> obj_cat -> count
        
    def to_dict(self) -> Dict[str, Any]:
        """Convert stats to dictionary for JSON serialization"""
        return {
            'source': self.source_name,
            'summary': {
                'total_nodes': self.total_nodes,
                'total_edges': self.total_edges,
                'unique_node_categories': len(self.node_categories),
                'unique_edge_predicates': len(self.edge_predicates),
                'unique_category_pairs': len(self.category_pairs)
            },
            'node_categories': dict(self.node_categories.most_common()),
            'edge_predicates': dict(self.edge_predicates.most_common()),
            'category_pairs': {
                f"{subj}--{obj}": count 
                for (subj, obj), count in self.category_pairs.most_common()
            },
            'predicate_category_distribution': {
                predicate: {
                    f"{subj}--{obj}": count
                    for (pred, subj, obj), count in self.predicate_category_pairs.items()
                    if pred == predicate
                }
                for predicate in self.edge_predicates.keys()
            }
        }


def generate_metagraph_streaming(nodes_file: Path, edges_file: Path, source_name: str = None) -> MetagraphStats:
    """Generate metagraph statistics from JSONL files using streaming"""
    
    if source_name is None:
        source_name = nodes_file.parent.name
    
    logging.info(f"Generating metagraph for {source_name}")
    
    stats = MetagraphStats(source_name)
    
    # Phase 1: Analyze nodes and build category mapping
    node_categories = {}  # node_id -> category
    
    logging.info("Analyzing nodes...")
    for node in stream_nodes_from_jsonl(nodes_file):
        node_id = node['id']
        category = normalize_category(node.get('category', 'biolink:NamedThing'))
        
        node_categories[node_id] = category
        stats.node_categories[category] += 1
        stats.total_nodes += 1
        
        if stats.total_nodes % 50000 == 0:
            logging.info(f"Processed {stats.total_nodes} nodes")
    
    logging.info(f"Found {len(stats.node_categories)} unique node categories")
    
    # Phase 2: Analyze edges
    logging.info("Analyzing edges...")
    for edge in stream_edges_from_jsonl(edges_file):
        subject_id = edge['subject']
        object_id = edge['object']
        predicate = edge.get('predicate', 'biolink:related_to')
        
        # Get categories (default to NamedThing if not found)
        subject_category = node_categories.get(subject_id, 'biolink:NamedThing')
        object_category = node_categories.get(object_id, 'biolink:NamedThing')
        
        # Update statistics
        stats.edge_predicates[predicate] += 1
        stats.category_pairs[(subject_category, object_category)] += 1
        stats.predicate_category_pairs[(predicate, subject_category, object_category)] += 1
        
        # Update node degrees
        stats.node_degrees[subject_id] += 1
        stats.node_degrees[object_id] += 1
        
        # Update category connectivity
        stats.category_connectivity[subject_category][object_category] += 1
        
        stats.total_edges += 1
        
        if stats.total_edges % 100000 == 0:
            logging.info(f"Processed {stats.total_edges} edges")
    
    logging.info(f"Metagraph analysis complete: {stats.total_nodes} nodes, {stats.total_edges} edges")
    return stats


def save_metagraph(stats: MetagraphStats, output_file: Path):
    """Save metagraph statistics to JSON file"""
    logging.info(f"Saving metagraph to {output_file}")
    
    with open(output_file, 'w') as f:
        json.dump(stats.to_dict(), f, indent=2)
    
    logging.info(f"Metagraph saved: {stats.total_nodes} nodes, {stats.total_edges} edges")


def generate_metagraph_for_source(nodes_file: Path, edges_file: Path, output_dir: Path, source_name: str = None):
    """Generate and save metagraph for a single source"""
    if source_name is None:
        source_name = nodes_file.parent.name
    
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{source_name}_metagraph.json"
    
    stats = generate_metagraph_streaming(nodes_file, edges_file, source_name)
    save_metagraph(stats, output_file)
    
    return output_file


def compare_metagraphs(metagraph_files: list, output_file: Path):
    """Compare multiple metagraphs and generate comparison report"""
    logging.info(f"Comparing {len(metagraph_files)} metagraphs")
    
    metagraphs = []
    for file_path in metagraph_files:
        with open(file_path, 'r') as f:
            metagraphs.append(json.load(f))
    
    comparison = {
        'sources_compared': [mg['source'] for mg in metagraphs],
        'summary_comparison': {},
        'category_overlap': {},
        'predicate_overlap': {},
        'unique_to_source': {}
    }
    
    # Summary comparison
    for mg in metagraphs:
        source = mg['source']
        comparison['summary_comparison'][source] = mg['summary']
    
    # Find overlaps and unique elements
    all_categories = set()
    all_predicates = set()
    
    for mg in metagraphs:
        categories = set(mg['node_categories'].keys())
        predicates = set(mg['edge_predicates'].keys())
        
        all_categories.update(categories)
        all_predicates.update(predicates)
        
        source = mg['source']
        comparison['unique_to_source'][source] = {
            'categories': list(categories),
            'predicates': list(predicates)
        }
    
    # Calculate overlaps
    comparison['category_overlap'] = {
        'total_unique_categories': len(all_categories),
        'categories': list(all_categories)
    }
    
    comparison['predicate_overlap'] = {
        'total_unique_predicates': len(all_predicates), 
        'predicates': list(all_predicates)
    }
    
    # Save comparison
    with open(output_file, 'w') as f:
        json.dump(comparison, f, indent=2)
    
    logging.info(f"Metagraph comparison saved to {output_file}")


def generate_metagraph_summary(stats: MetagraphStats) -> str:
    """Generate human-readable summary of metagraph"""
    summary_lines = [
        f"=== Metagraph Summary: {stats.source_name} ===",
        f"Total Nodes: {stats.total_nodes:,}",
        f"Total Edges: {stats.total_edges:,}",
        f"Unique Node Categories: {len(stats.node_categories)}",
        f"Unique Edge Predicates: {len(stats.edge_predicates)}",
        "",
        "Top Node Categories:",
    ]
    
    for category, count in stats.node_categories.most_common(10):
        percentage = (count / stats.total_nodes) * 100
        summary_lines.append(f"  {category}: {count:,} ({percentage:.1f}%)")
    
    summary_lines.extend([
        "",
        "Top Edge Predicates:",
    ])
    
    for predicate, count in stats.edge_predicates.most_common(10):
        percentage = (count / stats.total_edges) * 100
        summary_lines.append(f"  {predicate}: {count:,} ({percentage:.1f}%)")
    
    summary_lines.extend([
        "",
        "Top Category Pairs:",
    ])
    
    for (subj_cat, obj_cat), count in stats.category_pairs.most_common(10):
        percentage = (count / stats.total_edges) * 100
        summary_lines.append(f"  {subj_cat} -> {obj_cat}: {count:,} ({percentage:.1f}%)")
    
    return "\n".join(summary_lines)


def normalize_category(category: Any) -> str:
    """Normalize category to consistent format"""
    if isinstance(category, list):
        # Use first category if multiple
        category = category[0] if category else 'biolink:NamedThing'
    
    if not isinstance(category, str):
        return 'biolink:NamedThing'
    
    # Ensure biolink prefix
    if not category.startswith('biolink:'):
        category = f'biolink:{category}'
    
    return category


def create_cytoscape_metagraph(stats: MetagraphStats, output_file: Path, min_edge_count: int = 1):
    """Create Cytoscape-compatible metagraph visualization file"""
    logging.info(f"Creating Cytoscape metagraph for {stats.source_name}")
    
    # Create nodes (categories)
    nodes = []
    for category, count in stats.node_categories.items():
        nodes.append({
            'data': {
                'id': category,
                'label': category.replace('biolink:', ''),
                'node_count': count,
                'size': min(100, max(10, count // 1000))  # Scale node size
            }
        })
    
    # Create edges (category relationships)
    edges = []
    edge_id = 0
    for (source_cat, target_cat), count in stats.category_pairs.items():
        if count >= min_edge_count:
            edges.append({
                'data': {
                    'id': f'edge_{edge_id}',
                    'source': source_cat,
                    'target': target_cat,
                    'edge_count': count,
                    'weight': min(10, max(1, count // 1000))  # Scale edge weight
                }
            })
            edge_id += 1
    
    cytoscape_data = {
        'elements': {
            'nodes': nodes,
            'edges': edges
        },
        'metadata': {
            'source': stats.source_name,
            'total_nodes': stats.total_nodes,
            'total_edges': stats.total_edges,
            'min_edge_count_filter': min_edge_count
        }
    }
    
    with open(output_file, 'w') as f:
        json.dump(cytoscape_data, f, indent=2)
    
    logging.info(f"Cytoscape metagraph saved: {len(nodes)} category nodes, {len(edges)} category edges")