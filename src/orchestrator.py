"""
Main orchestration functions for PhenomeKG build
Streaming approach to workflow management
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers.kg2 import harmonize_kg2
from .harmonizers.spoke import harmonize_spoke
from .integration.entity_resolution import integrate_sources_streaming
from .post_processing.arango_export import prepare_for_arango_streaming
from .post_processing.biomapper_export import export_for_biomapper_streaming


def run_kg_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building PhenomeKG"""
    logging.info("Starting PhenomeKG build...")

    # Phase 1: Harmonize all sources to Biolink schema
    harmonized_sources = harmonize_all_sources(config['sources'], config.get('metagraph', {}))

    # Phase 2: Integrate into unified KG with entity resolution
    integration_config = config['integration'].copy()
    integration_config['metagraph_config'] = {
        'generate_summaries': config.get('metagraph', {}).get('generate_summaries', True),
        'generate_cytoscape': config.get('metagraph', {}).get('generate_cytoscape', True),
        'generate_html_viewer': config.get('metagraph', {}).get('generate_cytoscape', True),
        'cytoscape_thresholds': [1, 5, 10, 25]
    }
    
    unified_nodes, unified_edges = integrate_sources_streaming(
        harmonized_sources, 
        Path(config['integration']['output_directory']),
        integration_config
    )

    # Phase 3: Post-processing steps
    if 'post_processing' in config:
        post_process_unified_kg(unified_nodes, unified_edges, config['post_processing'])

    logging.info(f"Build complete: {unified_nodes}, {unified_edges}")
    return unified_nodes, unified_edges


def harmonize_all_sources(sources_config: dict, metagraph_config: dict) -> Dict[str, Dict[str, Path]]:
    """Harmonize each source that needs it"""
    harmonized_sources = {}

    for source_name, source_config in sources_config.items():
        if needs_harmonization(source_name, source_config):
            nodes_path, edges_path = harmonize_source(source_name, source_config, metagraph_config)
        else:
            logging.info(f"Skipping {source_name} - already harmonized")
            nodes_path = Path(source_config['harmonized_output']['nodes'])
            edges_path = Path(source_config['harmonized_output']['edges'])

        harmonized_sources[source_name] = {
            'nodes': nodes_path,
            'edges': edges_path
        }

    return harmonized_sources


def harmonize_source(source_name: str, config: dict, metagraph_config: dict) -> tuple[Path, Path]:
    """Harmonize a single source to Biolink schema"""
    logging.info(f"Harmonizing {source_name}...")

    # Get output paths
    nodes_output = Path(config['harmonized_output']['nodes'])
    edges_output = Path(config['harmonized_output']['edges'])

    # Create output directory if it doesn't exist
    nodes_output.parent.mkdir(parents=True, exist_ok=True)
    edges_output.parent.mkdir(parents=True, exist_ok=True)

    # Prepare harmonization rules
    harmonization_rules = config.get('harmonization_rules', {})
    harmonization_rules['generate_metagraph'] = metagraph_config.get('generate_for_sources', True)
    harmonization_rules['metagraph_config'] = {
        'generate_summaries': metagraph_config.get('generate_summaries', True),
        'generate_cytoscape': metagraph_config.get('generate_cytoscape', True),
        'generate_html_viewer': metagraph_config.get('generate_cytoscape', True),
        'cytoscape_thresholds': [1, 5, 10]
    }

    # Run source-specific harmonizer
    if source_name == 'kg2':
        harmonize_kg2(
            Path(config['nodes_input']),
            Path(config['edges_input']),
            nodes_output,
            edges_output,
            harmonization_rules
        )
    elif source_name == 'spoke':
        harmonize_spoke(
            Path(config['input_file']),
            nodes_output,
            edges_output,
            harmonization_rules
        )
    else:
        raise ValueError(f"Unknown source type: {source_name}")

    logging.info(f"Harmonized {source_name} -> {nodes_output}, {edges_output}")
    return nodes_output, edges_output


def needs_harmonization(source_name: str, config: dict) -> bool:
    """Check if harmonization step needs to run based on file timestamps"""
    # Get input paths based on source type
    if config.get('source_type') == 'separate_files':
        input_paths = [Path(config['nodes_input']), Path(config['edges_input'])]
    else:
        input_paths = [Path(config['input_file'])]
    
    # Get output paths
    nodes_output = Path(config['harmonized_output']['nodes'])
    edges_output = Path(config['harmonized_output']['edges'])
    output_paths = [nodes_output, edges_output]

    # Check if any output files are missing
    for output_path in output_paths:
        if not output_path.exists():
            logging.info(f"{source_name} needs harmonization: {output_path.name} doesn't exist")
            return True

    # Check if any input files are missing
    for input_path in input_paths:
        if not input_path.exists():
            logging.warning(f"{source_name} input file doesn't exist: {input_path}")
            return False

    # Check timestamps - need to harmonize if any input is newer than any output
    latest_input_mtime = max(path.stat().st_mtime for path in input_paths)
    earliest_output_mtime = min(path.stat().st_mtime for path in output_paths)

    needs_update = latest_input_mtime > earliest_output_mtime
    if needs_update:
        logging.info(f"{source_name} needs harmonization: input is newer than output")

    return needs_update


def post_process_unified_kg(unified_nodes: Path, unified_edges: Path, config: dict):
    """Run all post-processing steps on the unified KG"""
    logging.info("Starting post-processing...")

    # Step 1: Prepare ArangoDB version
    if 'arango_export' in config and config['arango_export'].get('enabled', True):
        arango_config = config['arango_export']
        output_dir = Path(arango_config['output_directory'])
        
        if needs_arango_export(unified_nodes, unified_edges, arango_config):
            logging.info("Preparing ArangoDB export...")
            prepare_for_arango_streaming(unified_nodes, unified_edges, output_dir, arango_config)
        else:
            logging.info("Skipping ArangoDB export - up to date")

    # Step 2: Export for biomapper
    if 'biomapper_export' in config and config['biomapper_export'].get('enabled', True):
        biomapper_config = config['biomapper_export']
        output_dir = Path(biomapper_config['output_directory'])
        
        if needs_biomapper_export(unified_nodes, unified_edges, biomapper_config):
            logging.info("Exporting for biomapper...")
            export_for_biomapper_streaming(unified_nodes, unified_edges, output_dir, biomapper_config)
        else:
            logging.info("Skipping biomapper export - up to date")

    logging.info("Post-processing complete")


def needs_arango_export(unified_nodes: Path, unified_edges: Path, config: dict) -> bool:
    """Check if ArangoDB export needs to run"""
    output_dir = Path(config['output_directory'])
    expected_files = config['output_files']
    
    nodes_output = output_dir / expected_files['nodes']
    edges_output = output_dir / expected_files['edges']

    # Check if output files exist
    if not nodes_output.exists() or not edges_output.exists():
        return True

    # Check timestamps
    unified_mtime = max(unified_nodes.stat().st_mtime, unified_edges.stat().st_mtime)
    output_mtime = min(nodes_output.stat().st_mtime, edges_output.stat().st_mtime)

    return unified_mtime > output_mtime


def needs_biomapper_export(unified_nodes: Path, unified_edges: Path, config: dict) -> bool:
    """Check if biomapper export needs to run"""
    output_dir = Path(config['output_directory'])

    if not output_dir.exists():
        return True

    # Check if any of the expected output files are missing or outdated
    unified_mtime = max(unified_nodes.stat().st_mtime, unified_edges.stat().st_mtime)

    # If output directory is empty, need to export
    output_files = list(output_dir.iterdir())
    if not output_files:
        return True

    # Check if any output files are older than unified files
    for file_path in output_files:
        if file_path.is_file() and file_path.stat().st_mtime < unified_mtime:
            return True

    return False