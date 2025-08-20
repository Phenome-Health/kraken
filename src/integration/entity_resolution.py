"""
Entity resolution and graph integration functions
"""

from pathlib import Path
from typing import Dict, Optional, List, Iterator
import logging
import json

from ..utils.kg_io import (
    stream_nodes_from_jsonl, 
    stream_edges_from_jsonl,
    load_equivalency_mappings,
    save_nodes_to_jsonl,
    save_edges_to_jsonl
)
from ..utils.metagraph import generate_metagraph_for_source, compare_metagraphs


def integrate_sources_streaming(harmonized_sources: Dict[str, Dict[str, Path]], output_dir: Path, config: dict):
    """Merge harmonized sources using streaming approach"""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output files
    unified_nodes = output_dir / "unified_nodes.jsonl"
    unified_edges = output_dir / "unified_edges.jsonl"
    
    logging.info("Starting streaming source integration...")
    
    # Phase 1: Build equivalency mappings from primary source
    primary_source = config.get('primary_source', 'rtx-kg2')
    primary_nodes_file = harmonized_sources[primary_source]['nodes']
    
    logging.info(f"Loading equivalency mappings from {primary_source}")
    equivalency_index = load_equivalency_mappings(primary_nodes_file)
    canonical_nodes = {}  # canonical_id -> merged_node_data
    
    # Phase 2: Process all nodes, merging as we go
    processed_nodes = set()
    
    for source_name, source_files in harmonized_sources.items():
        logging.info(f"Processing nodes from {source_name}")
        nodes_file = source_files['nodes']
        
        for node in stream_nodes_from_jsonl(nodes_file):
            node_id = node['id']
            
            # Find canonical ID for this node
            canonical_id = find_canonical_id(node_id, node.get('equivalent_identifiers', []), equivalency_index)
            
            if canonical_id in processed_nodes:
                # Merge with existing canonical node
                canonical_nodes[canonical_id] = merge_node_data(
                    canonical_nodes[canonical_id], node, source_name, config
                )
            else:
                # First time seeing this canonical entity
                canonical_nodes[canonical_id] = node.copy()
                canonical_nodes[canonical_id]['id'] = canonical_id
                canonical_nodes[canonical_id]['sources'] = [source_name]
                processed_nodes.add(canonical_id)
    
    logging.info(f"Processed {len(canonical_nodes)} unique nodes")
    
    # Save unified nodes
    save_nodes_to_jsonl(canonical_nodes.values(), unified_nodes)
    
    # Phase 3: Process all edges with node ID resolution
    edge_dedup = set()  # (canonical_subject, canonical_object, predicate) for deduplication
    
    def process_edges():
        for source_name, source_files in harmonized_sources.items():
            logging.info(f"Processing edges from {source_name}")
            edges_file = source_files['edges']
            
            for edge in stream_edges_from_jsonl(edges_file):
                # Resolve subject and object to canonical IDs
                canonical_subject = find_canonical_id(
                    edge['subject'], [], equivalency_index
                )
                canonical_object = find_canonical_id(
                    edge['object'], [], equivalency_index
                )
                
                # Create edge key for deduplication
                predicate = edge.get('predicate', 'biolink:related_to')
                edge_key = (canonical_subject, canonical_object, predicate)
                
                if edge_key not in edge_dedup:
                    edge_dedup.add(edge_key)
                    
                    # Update edge with canonical IDs
                    unified_edge = edge.copy()
                    unified_edge['subject'] = canonical_subject
                    unified_edge['object'] = canonical_object
                    unified_edge['sources'] = [source_name]
                    
                    yield unified_edge
    
    # Save unified edges
    save_edges_to_jsonl(process_edges(), unified_edges)
    
    logging.info(f"Integration complete! Unified KG saved to {output_dir}")
    
    # Generate metagraph for unified result
    if config.get('generate_metagraph', True):
        # Store unified metagraphs in artifacts/metagraphs/unified/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "unified"
        
        metagraph_config = config.get('metagraph_config', {
            'generate_summaries': True,
            'generate_cytoscape': True,
            'generate_html_viewer': True,
            'cytoscape_thresholds': [1, 5, 10, 25]  # Additional threshold for unified graph
        })
        
        unified_metagraph_files = generate_metagraph_for_source(
            unified_nodes, unified_edges, metagraph_dir, "unified", metagraph_config
        )
        logging.info("Unified metagraph generated")
        
        # Compare with source metagraphs if they exist
        source_metagraphs = []
        for source_name in harmonized_sources.keys():
            source_metagraph = artifacts_root / "metagraphs" / "harmonized" / source_name / f"{source_name}_metagraph.json"
            if source_metagraph.exists():
                source_metagraphs.append(source_metagraph)
        
        if source_metagraphs:
            # Find the main JSON file from unified metagraph
            unified_json = next((f for f in unified_metagraph_files if f.name.endswith('_metagraph.json')), None)
            if unified_json:
                source_metagraphs.append(unified_json)
                comparison_file = metagraph_dir / "metagraph_comparison.json"
                compare_metagraphs(source_metagraphs, comparison_file)
                logging.info("Metagraph comparison generated")
    
    return unified_nodes, unified_edges


def find_canonical_id(node_id: str, equiv_ids: List[str], equivalency_index: Dict[str, set]) -> str:
    """Find canonical ID for this node using equivalency mappings"""
    # Check if this node or any equivalent ID is in our index
    all_ids = [node_id] + equiv_ids
    
    for id_val in all_ids:
        if id_val in equivalency_index:
            # Return the first (canonical) ID from the equivalency set
            return min(equivalency_index[id_val])
    
    # Not found in index, use the original ID as canonical
    return node_id


def merge_source_with_resolution(base_kg, source_kg, node_index, source_name, config):
    """Merge source KG into base with entity resolution"""
    namespace_priority = config.get('entity_resolution', {}).get('namespace_priority', [])

    nodes_merged = 0
    nodes_added = 0

    # Process nodes
    for source_node_id, source_node_data in source_kg.nodes(data=True):
        source_equiv_ids = source_node_data.get('equivalent_ids', [])

        # Check if this node matches any existing node
        canonical_id = find_matching_node(source_node_id, source_equiv_ids, node_index)

        if canonical_id:
            # Merge with existing node
            existing_node_data = base_kg.nodes[canonical_id]
            merged_node_data = merge_node_data(
                existing_node_data, source_node_data, source_name, namespace_priority
            )
            base_kg.nodes[canonical_id].update(merged_node_data)

            # Update index with any new equivalent IDs
            for equiv_id in source_equiv_ids:
                if equiv_id not in node_index:
                    node_index[equiv_id] = canonical_id

            nodes_merged += 1
        else:
            # Add as new node
            base_kg.add_node(source_node_id, **source_node_data)

            # Update index
            node_index[source_node_id] = source_node_id
            for equiv_id in source_equiv_ids:
                node_index[equiv_id] = source_node_id

            nodes_added += 1

    # Process edges (with node ID resolution)
    edges_merged = 0
    edges_added = 0

    for source_edge in source_kg.edges(data=True):
        source_subject, source_object, edge_data = source_edge

        # Resolve subject and object to canonical IDs
        canonical_subject = node_index.get(source_subject, source_subject)
        canonical_object = node_index.get(source_object, source_object)

        # Check if edge already exists
        if base_kg.has_edge(canonical_subject, canonical_object):
            # Merge edge data
            existing_edge_data = base_kg.edges[canonical_subject, canonical_object]
            merged_edge_data = merge_edge_data(existing_edge_data, edge_data, source_name)
            base_kg.edges[canonical_subject, canonical_object].update(merged_edge_data)
            edges_merged += 1
        else:
            # Add new edge
            base_kg.add_edge(canonical_subject, canonical_object, **edge_data)
            edges_added += 1

    logging.info(f"  Nodes: {nodes_merged} merged, {nodes_added} added")
    logging.info(f"  Edges: {edges_merged} merged, {edges_added} added")

    return base_kg, node_index


def find_matching_node(node_id: str, equiv_ids: List[str], node_index: Dict[str, str]) -> Optional[str]:
    """Find if this node matches any existing node via equivalent IDs"""
    all_ids = [node_id] + equiv_ids

    for id_val in all_ids:
        if id_val in node_index:
            return node_index[id_val]

    return None


def merge_node_data(existing_node: dict, new_node: dict, source_name: str, config: dict) -> dict:
    """Merge data from new node into existing node"""
    merged = existing_node.copy()

    # Merge equivalent_ids
    existing_equivs = set(merged.get('equivalent_identifiers', []))
    new_equivs = set(new_node.get('equivalent_identifiers', []))
    merged['equivalent_identifiers'] = list(existing_equivs | new_equivs)

    # Add source provenance  
    existing_sources = merged.get('sources', [])
    if source_name not in existing_sources:
        existing_sources.append(source_name)
    merged['sources'] = existing_sources

    # Merge other properties (simple first-wins strategy for now)
    for key, value in new_node.items():
        if key == 'id':
            # Keep existing canonical ID
            continue
        elif key not in merged or merged[key] is None:
            # Add new property
            merged[key] = value
        # For conflicts, keep existing value (could be made more sophisticated)

    return merged


def merge_edge_data(existing_edge: dict, new_edge: dict, source_name: str) -> dict:
    """Merge edge data, combining sources and handling conflicts"""
    merged = existing_edge.copy()

    # Add source provenance
    sources = set(merged.get('sources', [existing_edge.get('source', '')]))
    sources.add(source_name)
    merged['sources'] = list(sources)

    # For edges, we mostly just want to track multiple sources
    # Could add more sophisticated conflict resolution here if needed

    return merged


def resolve_property_conflict(existing_value, new_value, property_name: str, namespace_priority: List[str]):
    """Resolve conflicts between property values from different sources"""

    # If values are the same, no conflict
    if existing_value == new_value:
        return existing_value

    # For now, just keep the existing value (first source wins)
    # Could implement more sophisticated resolution based on source priority
    logging.debug(f"Property conflict for {property_name}: keeping '{existing_value}' over '{new_value}'")
    return existing_value


def choose_canonical_id(equivalent_ids: set, namespace_priority: List[str] = None) -> str:
    """Choose the canonical ID from a set of equivalent IDs"""
    if not namespace_priority:
        namespace_priority = [
            'CHEMBL.COMPOUND',
            'PUBCHEM.COMPOUND',
            'DRUGBANK',
            'MESH',
            'UNIPROT',
            'ENSEMBL',
            'HGNC'
        ]

    # Try each namespace in priority order
    for namespace in namespace_priority:
        candidates = [id_val for id_val in equivalent_ids if id_val.startswith(namespace)]
        if candidates:
            return sorted(candidates)[0]  # Take first alphabetically if multiple

    # Fallback: shortest ID or lexicographically first
    return min(equivalent_ids, key=lambda x: (len(x), x))
