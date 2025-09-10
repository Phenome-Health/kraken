"""
Main orchestration functions for PhenomeKG build
Streaming approach to workflow management
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers.kg2 import harmonize_kg2
from .harmonizers.spoke import harmonize_spoke
from .harmonizers.umls import harmonize_umls
from .integration.entity_resolution import integrate_sources
from .post_processing.arango_export import export_for_arango
from .post_processing.biomapper_export import export_for_biomapper
from .utils.metagraph import generate_metagraph_for_source, compare_metagraphs


def run_kg_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building PhenomeKG"""
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
        integrate_sources(harmonized_source_paths, unified_dir_path, unified_nodes_path, unified_edges_path, config['integration'].copy())

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
    if source_name == 'kg2':
        harmonize_kg2(
            Path(config['nodes_input']),
            Path(config['edges_input']),
            nodes_output,
            edges_output,
            biolink_version,
            build_metagraph
        )
    elif source_name == 'spoke':
        harmonize_spoke(
            Path(config['input_file']),
            nodes_output,
            edges_output,
            biolink_version,
            build_metagraph
        )
    elif source_name == 'umls':
        harmonize_umls(
            Path(config['input_file']),
            nodes_output,
            edges_output,
            biolink_version,
            build_metagraph
        )
    else:
        raise ValueError(f"Unknown source type: {source_name}")

    logging.info(f"Harmonized {source_name} -> {nodes_output}, {edges_output}")


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
    arango_output_dir = Path(config['arango_export']['output_directory'])
    arango_nodes_path = arango_output_dir / f"kraken_{kraken_version}_nodes_arango.jsonl"
    arango_edges_path = arango_output_dir / f"kraken_{kraken_version}_edges_arango.jsonl"

    # Step 1: Prepare ArangoDB version
    if config['arango_export'].get('enabled', True):
        logging.info("Preparing ArangoDB export...")
        export_for_arango(unified_nodes_path, unified_edges_path, arango_output_dir,
                          arango_nodes_path, arango_edges_path, biolink_version, kraken_version)

    # Step 2: Export for biomapper (off of Arango export files)
    if 'biomapper_export' in config and config['biomapper_export'].get('enabled', True):
        biomapper_config = config['biomapper_export']
        output_dir = Path(biomapper_config['output_directory'])
        
        logging.info("Exporting for biomapper...")
        export_for_biomapper(arango_nodes_path, output_dir, kraken_version)

    logging.info("Post-processing complete")
