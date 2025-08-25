"""
Main orchestration functions for PhenomeKG build
Streaming approach to workflow management
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers.kg2 import harmonize_kg2
from .harmonizers.spoke import harmonize_spoke
from .integration.entity_resolution import integrate_sources
from .post_processing.arango_export import prepare_for_arango
from .post_processing.biomapper_export import export_for_biomapper
from .utils.metagraph import generate_metagraph_for_source, compare_metagraphs


def run_kg_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building PhenomeKG"""
    biolink_version = config['biolink_version']
    unified_dir_path = Path(config['integration']['output_directory'])
    unified_nodes_path = unified_dir_path / config['integration']['unified_output']['nodes']
    unified_edges_path = unified_dir_path / config['integration']['unified_output']['edges']
    harmonized_source_paths = {source_name: {'nodes': source_config['harmonized_output']['nodes'],
                                             'edges': source_config['harmonized_output']['edges']} 
                                for source_name, source_config in config['sources'].items()}

    # Phase 1: Harmonize all sources to Biolink semantic layer/schema
    if config['steps'].get('harmonize'):
        harmonize_all_sources(config['sources'], biolink_version)

    # Phase 2: Integrate into unified KG with entity resolution
    if config['steps'].get('integrate'):
        integrate_sources(harmonized_source_paths, unified_dir_path, config['integration'].copy())

    # Phase 3: Generate metagraph for unified result
    if config['steps'].get('metagraph'):
        generate_unified_metagraph(unified_nodes_path, unified_edges_path, harmonized_source_paths)

    # Phase 4: Post-processing steps
    if config['steps'].get('postprocess'):
        post_process_unified_kg(unified_nodes_path, unified_edges_path, config['post_processing'], biolink_version)

    logging.info(f"Build complete: {unified_nodes_path}, {unified_edges_path}")
    return unified_nodes_path, unified_edges_path


def harmonize_all_sources(sources_config: dict, biolink_version: str):
    """Harmonize each source that needs it"""
    for source_name, source_config in sources_config.items():
        # NOTE: for now, always re-harmonize with every build
        nodes_path, edges_path = harmonize_source(source_name, source_config, biolink_version)


def harmonize_source(source_name: str, config: dict, biolink_version: str) -> tuple[Path, Path]:
    """Harmonize a single source to Biolink schema"""
    logging.info(f"Harmonizing {source_name}...")

    # Get output paths
    nodes_output = Path(config['harmonized_output']['nodes'])
    edges_output = Path(config['harmonized_output']['edges'])

    # Create output directory if it doesn't exist
    nodes_output.parent.mkdir(parents=True, exist_ok=True)
    edges_output.parent.mkdir(parents=True, exist_ok=True)

    # Run source-specific harmonizer
    if source_name == 'kg2':
        harmonize_kg2(
            Path(config['nodes_input']),
            Path(config['edges_input']),
            nodes_output,
            edges_output,
            biolink_version
        )
    elif source_name == 'spoke':
        harmonize_spoke(
            Path(config['input_file']),
            nodes_output,
            edges_output,
            biolink_version
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


def generate_unified_metagraph(unified_nodes_path: str, unified_edges_path: str, harmonized_source_paths: dict):
    # Store unified metagraphs in artifacts/metagraphs/unified/
    artifacts_root = Path("artifacts")
    metagraph_dir = artifacts_root / "metagraphs" / "unified"
    
    unified_metagraph_files = generate_metagraph_for_source(unified_nodes_path, unified_edges_path, metagraph_dir, "unified")
    logging.info("Unified metagraph generated")
    
    # Compare with source metagraphs if they exist
    source_metagraphs = []
    for source_name in harmonized_source_paths.keys():
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


def post_process_unified_kg(unified_nodes_path: Path, unified_edges_path: Path, config: dict, biolink_version: str):
    """Run all post-processing steps on the unified KG"""
    logging.info("Starting post-processing...")

    # Step 1: Prepare ArangoDB version
    if 'arango_export' in config and config['arango_export'].get('enabled', True):
        arango_config = config['arango_export']
        output_dir = Path(arango_config['output_directory'])
        
        if needs_arango_export(unified_nodes_path, unified_edges_path, arango_config):
            logging.info("Preparing ArangoDB export...")
            arango_unified_nodes_path, arango_unified_edges_path = prepare_for_arango(unified_nodes_path, unified_edges_path, output_dir, arango_config, biolink_version)
        else:
            logging.info("Skipping ArangoDB export - up to date")

    # Step 2: Export for biomapper (off of Arango export files)
    if 'biomapper_export' in config and config['biomapper_export'].get('enabled', True):
        biomapper_config = config['biomapper_export']
        output_dir = Path(biomapper_config['output_directory'])
        
        if needs_biomapper_export(arango_unified_nodes_path, biomapper_config):
            logging.info("Exporting for biomapper...")
            export_for_biomapper(arango_unified_nodes_path, output_dir)
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


def needs_biomapper_export(unified_nodes: Path, config: dict) -> bool:
    """Check if biomapper export needs to run"""
    output_dir = Path(config['output_directory'])

    if not output_dir.exists():
        return True

    # Check if any of the expected output files are missing or outdated
    unified_mtime = unified_nodes.stat().st_mtime

    # If output directory is empty, need to export
    output_files = list(output_dir.iterdir())
    if not output_files:
        return True

    # Check if any output files are older than unified files
    for file_path in output_files:
        if file_path.is_file() and file_path.stat().st_mtime < unified_mtime:
            return True

    return False