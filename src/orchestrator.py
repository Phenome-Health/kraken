"""
Main orchestration functions for KRAKEN build
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers.clingen import harmonize_clingen
from .harmonizers.kg2 import harmonize_kg2
from .harmonizers.lipidmaps import harmonize_lipidmaps
from .harmonizers.refmet import harmonize_refmet
from .harmonizers.robokop import harmonize_robokop
from .harmonizers.spoke import harmonize_spoke
from .harmonizers.umls import harmonize_umls
from .integration.entity_resolution import integrate_sources
from .utils.kg_io import unzip_files, zip_files, get_harmonized_file_paths, PROJECT_ROOT
from .utils.metagraph import generate_metagraph_for_source, compare_metagraphs
from .utils.general import to_list
from .post_processing.test_file_generator import create_test_kg_files





def run_kg_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building the KRAKEN"""
    biolink_version = config['biolink_version']
    kraken_version = config['kraken_version']
    unified_dir_path = Path(config['integration']['output_directory'])
    unified_nodes_path = unified_dir_path / f"kraken_nodes_{kraken_version}.jsonl"
    unified_edges_path = unified_dir_path / f"kraken_edges_{kraken_version}.jsonl"
    user_specified_sources = to_list(config["options"].get("include_sources"))
    sources_to_use = user_specified_sources if user_specified_sources else list(config["sources"].keys())
    logging.info(f"Will include {len(sources_to_use)} sources: {sources_to_use}")

    create_metagraphs = True if config['options'].get('metagraph_creation') else False

    # Phase 1: Harmonize all sources to Biolink semantic layer/schema
    if config['steps'].get('harmonize'):
        logging.info(f"-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        harmonize_sources(config['sources'], biolink_version, build_metagraph=create_metagraphs)

    # Phase 2: Integrate into unified KG with entity resolution
    if config['steps'].get('integrate'):
        logging.info(f"-------------------------- INTEGRATING SOURCES -----------------------------------------------")
        integrate_sources(sources_to_use, unified_dir_path, unified_nodes_path, unified_edges_path, config)

        if create_metagraphs:
            logging.info(f"---------------------- GENERATING UNIFIED METAGRAPH ------------------------------------------")
            generate_unified_metagraph(unified_nodes_path, unified_edges_path, sources_to_use)

    # Phase 3: Post-processing steps
    if config['steps'].get('postprocess'):
        logging.info(f"------------------------------ POST-PROCESSING -----------------------------------------------")
        post_process_unified_kg(unified_nodes_path, unified_edges_path, config['post_processing'], biolink_version, kraken_version)

    logging.info(f"Build complete: {unified_nodes_path}, {unified_edges_path}")
    return unified_nodes_path, unified_edges_path


def harmonize_sources(sources_config: dict, biolink_version: str, build_metagraph: bool):
    """Harmonize each source that needs it"""
    for source_name, source_config in sources_config.items():
        # NOTE: for now, always re-harmonize with every build
        harmonize_source(source_name, source_config, biolink_version, build_metagraph)


def harmonize_source(source_name: str, config: dict, biolink_version: str, build_metagraph: bool):
    """Harmonize a single source to Biolink schema"""
    logging.info(f"Harmonizing {source_name}...")

    # Get output paths
    nodes_output, edges_output = get_harmonized_file_paths(source_name)

    # Create output directory if it doesn't exist
    nodes_output.parent.mkdir(parents=True, exist_ok=True)
    edges_output.parent.mkdir(parents=True, exist_ok=True)

    # Run source-specific harmonizer
    harmonizers = {
        'clingen': harmonize_clingen,
        'kg2': harmonize_kg2,
        'robokop': harmonize_robokop,
        'spoke': harmonize_spoke,
        'umls': harmonize_umls,
        'lipidmaps': harmonize_lipidmaps,
        'refmet': harmonize_refmet
    }
    harmonizer = harmonizers[source_name]
    possible_input_file_fields = ['input_file', 'input_nodes', 'input_edges']
    input_file_paths = [config.get(field) for field in possible_input_file_fields]

    unzip_files(input_file_paths)

    if config.get('input_file'):
        harmonizer(config['input_file'], nodes_output, edges_output, biolink_version)
    elif config.get('nodes_input') and config.get('edges_input'):
        harmonizer(config['nodes_input'], config['edges_input'], nodes_output, edges_output, biolink_version)
    else:
        raise ValueError(f"Unknown source type: {source_name}")

    logging.info(f"Harmonized {source_name} -> {nodes_output}, {edges_output}")

    zip_files(input_file_paths)

    if build_metagraph:
        # Generate metagraph for harmonized output, stored in artifacts/metagraphs/harmonized/<source_name>/
        artifacts_root = PROJECT_ROOT / "artifacts"
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / source_name
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, source_name)


def generate_unified_metagraph(unified_nodes_path: Path, unified_edges_path: Path, source_names: list[str]):
    # Store unified metagraphs in artifacts/metagraphs/unified/
    artifacts_root = PROJECT_ROOT / "artifacts"
    metagraph_dir = artifacts_root / "metagraphs" / "unified"
    
    unified_metagraph_files = generate_metagraph_for_source(unified_nodes_path, unified_edges_path, metagraph_dir, "unified")
    logging.info("Unified metagraph generated")
    
    # Compare with source metagraphs if they exist
    source_metagraphs = []
    for source_name in source_names:
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
