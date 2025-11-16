"""
Main orchestration functions for KRAKEN build
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers.acmg import harmonize_acmg
from .harmonizers.kg2 import harmonize_kg2
from .harmonizers.lipidmaps import harmonize_lipidmaps
from .harmonizers.refmet import harmonize_refmet
from .harmonizers.spoke import harmonize_spoke
from .harmonizers.umls import harmonize_umls
from .integration.entity_resolution import integrate_sources
from .utils.metagraph import generate_metagraph_for_source, compare_metagraphs
from .post_processing.test_file_generator import create_test_kg_files


def run_kg_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building the KRAKEN"""
    biolink_version = config['biolink_version']
    kraken_version = config['kraken_version']
    unified_dir_path = Path(config['integration']['output_directory'])
    unified_nodes_path = unified_dir_path / f"kraken_nodes_{kraken_version}.jsonl"
    unified_edges_path = unified_dir_path / f"kraken_edges_{kraken_version}.jsonl"
    harmonized_source_paths = {source_name: {'nodes': Path(source_config['harmonized_output']['nodes']),
                                             'edges': Path(source_config['harmonized_output']['edges'])}
                                for source_name, source_config in config['sources'].items()}


    # Phase 1: Harmonize all sources to Biolink semantic layer/schema
    if config['steps'].get('harmonize'):
        logging.info(f"-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        harmonize_all_sources(config['sources'], biolink_version, build_metagraph=config['steps']['metagraph'])

    # Phase 2: Integrate into unified KG with entity resolution
    if config['steps'].get('integrate'):
        logging.info(f"-------------------------- INTEGRATING SOURCES -----------------------------------------------")
        integrate_sources(harmonized_source_paths, unified_dir_path, unified_nodes_path, unified_edges_path, config)

    # Phase 3: Generate metagraph for unified result
    if config['steps'].get('metagraph'):
        logging.info(f"---------------------- GENERATING UNIFIED METAGRAPH ------------------------------------------")
        generate_unified_metagraph(unified_nodes_path, unified_edges_path, harmonized_source_paths)

    # Phase 4: Post-processing steps
    if config['steps'].get('postprocess'):
        logging.info(f"------------------------------ POST-PROCESSING -----------------------------------------------")
        post_process_unified_kg(unified_nodes_path, unified_edges_path, config['post_processing'], biolink_version, kraken_version)

    logging.info(f"Build complete: {unified_nodes_path}, {unified_edges_path}")
    return unified_nodes_path, unified_edges_path


def harmonize_all_sources(sources_config: dict, biolink_version: str, build_metagraph: bool):
    """Harmonize each source that needs it"""
    for source_name, source_config in sources_config.items():
        if source_config.get('enabled'):
            # NOTE: for now, always re-harmonize with every build
            harmonize_source(source_name, source_config, biolink_version, build_metagraph)


def harmonize_source(source_name: str, config: dict, biolink_version: str, build_metagraph: bool):
    """Harmonize a single source to Biolink schema"""
    logging.info(f"Harmonizing {source_name}...")

    # Get output paths
    nodes_output = Path(config['harmonized_output']['nodes'])
    edges_output = Path(config['harmonized_output']['edges'])

    # Create output directory if it doesn't exist
    nodes_output.parent.mkdir(parents=True, exist_ok=True)
    edges_output.parent.mkdir(parents=True, exist_ok=True)

    # Run source-specific harmonizer
    harmonizers = {
        'acmg': harmonize_acmg,
        'kg2': harmonize_kg2,
        'spoke': harmonize_spoke,
        'umls': harmonize_umls,
        'lipidmaps': harmonize_lipidmaps,
        'refmet': harmonize_refmet
    }
    harmonizer = harmonizers[source_name]
    if config.get('input_file'):
        harmonizer(Path(config['input_file']), nodes_output, edges_output, biolink_version)
    elif config.get('nodes_input') and config.get('edges_input'):
        harmonizer(Path(config['nodes_input']), Path(config['edges_input']), nodes_output, edges_output, biolink_version)
    else:
        raise ValueError(f"Unknown source type: {source_name}")

    logging.info(f"Harmonized {source_name} -> {nodes_output}, {edges_output}")

    if build_metagraph:
        # Generate metagraph for harmonized output, stored in artifacts/metagraphs/harmonized/<source_name>/
        artifacts_root = Path("artifacts")
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / source_name
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, source_name)


def generate_unified_metagraph(unified_nodes_path: Path, unified_edges_path: Path, harmonized_source_paths: dict):
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


def post_process_unified_kg(unified_nodes_path: Path, unified_edges_path: Path, config: dict, biolink_version: str, kraken_version: str):
    """Run all post-processing steps on the unified KG"""
    logging.info("Starting post-processing...")

    if config.get('test_export'):
        logging.info("Generating test files for this kraken build..")
        test_export_config = config['test_export']
        output_dir = Path(test_export_config['output_directory'])
        create_test_kg_files(unified_nodes_path, unified_edges_path, output_dir, num_edges=test_export_config['num_edges'])

    logging.info("Post-processing complete")
