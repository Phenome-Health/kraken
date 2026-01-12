"""
Main orchestration functions for KRAKEN build
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers.clingen import ClinGenHarmonizer
from .harmonizers.molepro import MoleProHarmonizer
from .harmonizers.microbiome_kg import MicrobiomeKGHarmonizer
from .harmonizers.multiomics_kg import MultiomicsKGHarmonizer
from .harmonizers.refmet import RefMetHarmonizer
from .harmonizers.robokop import RobokopHarmonizer
from .harmonizers.kg2 import KG2Harmonizer
from .harmonizers.lipidmaps import LipidMapsHarmonizer
from .harmonizers.spoke import SpokeHarmonizer
from .harmonizers.umls import UMLSHarmonizer

from .integration.entity_resolution import integrate_sources
from .utils.biolink_client import BiolinkClient
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

    # Figure out which sources to use based on the build config settings
    all_sources = set(config["sources"])
    include_sources = set(to_list(config["options"].get("include_sources")))
    exclude_sources = set(to_list(config["options"].get("exclude_sources")))
    if not include_sources.issubset(all_sources) or not exclude_sources.issubset(all_sources):
        raise ValueError(f"In build_config.yaml, source names specified in 'include_sources' and "
                         f"'exclude_sources' must exist under the 'sources' slot, "
                         f"which currently includes only these: {all_sources}")
    overlap = include_sources.intersection(exclude_sources)
    if overlap:
        raise ValueError(f"Sources cannot be both included and excluded in build_config.yaml: {overlap}")
    source_pool = include_sources if include_sources else all_sources
    sources_to_use = source_pool - exclude_sources
    logging.info(f"Will include {len(sources_to_use)} sources: {sources_to_use}")

    create_metagraphs = True if config['options'].get('metagraph_creation') else False
    zip_inputs_after = True if config['options'].get('zip_inputs_after') else False

    # Phase 1: Harmonize all sources to Biolink semantic layer/schema
    source_configs = {source: config["sources"][source] for source in sources_to_use}
    if config['steps'].get('harmonize'):
        logging.info(f"-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        harmonize_sources(source_configs, biolink_version, create_metagraphs, zip_inputs_after)

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

    return unified_nodes_path, unified_edges_path


def harmonize_sources(sources_config: dict, biolink_version: str, build_metagraph: bool, zip_inputs_after: bool):
    """Harmonize each source that needs it"""
    biolink_client = BiolinkClient(biolink_version)
    for source_name, source_config in sources_config.items():
        # NOTE: for now, always re-harmonize with every build
        harmonize_source(source_name, source_config, biolink_client, build_metagraph, zip_inputs_after)


def harmonize_source(source_name: str, config: dict, biolink_client: BiolinkClient, build_metagraph: bool, zip_inputs_after: bool):
    """Harmonize a single source to Biolink schema"""
    logging.info(f"Harmonizing {source_name}...")

    # Get output paths
    nodes_output, edges_output = get_harmonized_file_paths(source_name)

    # Create output directory if it doesn't exist
    nodes_output.parent.mkdir(parents=True, exist_ok=True)
    edges_output.parent.mkdir(parents=True, exist_ok=True)

    # Run source-specific harmonizer
    harmonizers = {
        'kg2': KG2Harmonizer,
        'robokop': RobokopHarmonizer,
        'molepro': MoleProHarmonizer,
        'microbiome-kg': MicrobiomeKGHarmonizer,
        'multiomics-kg': MultiomicsKGHarmonizer,
        'spoke': SpokeHarmonizer,
        'umls': UMLSHarmonizer,
        'lipidmaps': LipidMapsHarmonizer,
        'refmet': RefMetHarmonizer,
        'clingen': ClinGenHarmonizer,
    }

    # Instantiate our harmonizer
    harmonizer = harmonizers[source_name](biolink_client)

    # Unzip input files as needed
    possible_input_file_fields = ['input_file', 'nodes_input', 'edges_input']
    input_file_paths = [config.get(field) for field in possible_input_file_fields]
    unzip_files(input_file_paths)

    # Harmonize input files
    if config.get('input_file'):
        # Note: These are not compatible with the BaseHarmonizer - have separate harmonize() implementations
        harmonizer.harmonize(config['input_file'], nodes_output, edges_output)
    elif config.get('nodes_input') and config.get('edges_input'):
        # Note: These (with split nodes/edges files) use the BaseHarmonizer
        harmonizer.harmonize(config['nodes_input'], config['edges_input'], nodes_output, edges_output)
    else:
        raise ValueError(f"Unknown source type: {source_name}")

    if zip_inputs_after:
        # Zip input files back up
        zip_files(input_file_paths)

    if build_metagraph:
        # Generate metagraph for harmonized output, stored in artifacts/metagraphs/harmonized/<source_name>/
        artifacts_root = PROJECT_ROOT / "artifacts"
        metagraph_dir = artifacts_root / "metagraphs" / "harmonized" / source_name
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir, source_name)


def generate_unified_metagraph(unified_nodes_path: Path, unified_edges_path: Path, source_names: set[str]):
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
