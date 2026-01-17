"""
Main orchestration functions for KRAKEN build
"""

import logging
from pathlib import Path

from kraken.harmonizers.clingen import ClinGenHarmonizer
from kraken.harmonizers.kg2 import KG2Harmonizer
from kraken.harmonizers.lipidmaps import LipidMapsHarmonizer
from kraken.harmonizers.microbiome_kg import MicrobiomeKGHarmonizer
from kraken.harmonizers.molepro import MoleProHarmonizer
from kraken.harmonizers.multiomics_kg import MultiomicsKGHarmonizer
from kraken.harmonizers.refmet import RefMetHarmonizer
from kraken.harmonizers.robokop import RobokopHarmonizer
from kraken.harmonizers.spoke import SpokeHarmonizer
from kraken.harmonizers.umls import UMLSHarmonizer
from kraken.integration.entity_resolution import integrate_sources
from kraken.post_processing.test_file_generator import create_test_kg_files
from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.constants import PROJECT_ROOT
from kraken.utils.general import to_list
from kraken.utils.kg_io import get_harmonized_file_paths, unzip_files, zip_files
from kraken.utils.metagraph import generate_metagraph_for_source


def run_kg_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building the KRAKEN"""
    biolink_version = config["biolink_version"]
    kraken_version = config["kraken_version"]
    base_path = Path(config["base_path"]) if config.get("base_path") else PROJECT_ROOT

    harmonized_dir_path = base_path / Path(config["harmonization"]["output_directory"])
    metagraph_dir_path = base_path / Path(config["metagraph"]["output_directory"])
    integrated_dir_path = base_path / Path(config["integration"]["output_directory"])
    integrated_nodes_path = integrated_dir_path / f"kraken_nodes_{kraken_version}.jsonl"
    integrated_edges_path = integrated_dir_path / f"kraken_edges_{kraken_version}.jsonl"

    # Figure out which sources to use based on the build config settings
    all_sources = set(config["sources"])
    include_sources = set(to_list(config["options"].get("include_sources")))
    exclude_sources = set(to_list(config["options"].get("exclude_sources")))
    if not include_sources.issubset(all_sources) or not exclude_sources.issubset(all_sources):
        raise ValueError(
            f"In build_config.yaml, source names specified in 'include_sources' and "
            f"'exclude_sources' must exist under the 'sources' slot, "
            f"which currently includes only these: {all_sources}"
        )
    overlap = include_sources.intersection(exclude_sources)
    if overlap:
        raise ValueError(f"Sources cannot be both included and excluded in build_config.yaml: {overlap}")
    source_pool = include_sources if include_sources else all_sources
    sources_to_use = source_pool - exclude_sources
    logging.info(f"Will include {len(sources_to_use)} sources: {sources_to_use}")

    create_metagraphs = True if config["options"].get("metagraph_creation") else False
    zip_inputs_after = True if config["harmonization"].get("zip_inputs_after") else False

    # Phase 1: Harmonize all sources to Biolink semantic layer/schema
    source_configs = {source: config["sources"][source] for source in sources_to_use}
    if config["steps"].get("harmonize"):
        logging.info("-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        harmonize_sources(
            source_configs,
            biolink_version,
            base_path,
            harmonized_dir_path,
            create_metagraphs,
            metagraph_dir_path,
            zip_inputs_after,
        )

    # Phase 2: Integrate into unified KG with entity resolution
    if config["steps"].get("integrate"):
        logging.info("-------------------------- INTEGRATING SOURCES -----------------------------------------------")
        integrate_sources(
            sources_to_use,
            integrated_dir_path,
            integrated_nodes_path,
            integrated_edges_path,
            harmonized_dir_path,
            config,
        )

        if create_metagraphs:
            generate_metagraph_for_source(integrated_nodes_path, integrated_edges_path, metagraph_dir_path, "kraken")

    # Phase 3: Post-processing steps
    if config["steps"].get("postprocess"):
        logging.info("------------------------------ POST-PROCESSING -----------------------------------------------")
        post_process_integrated_kg(
            integrated_nodes_path,
            integrated_edges_path,
            base_path,
            config["post_processing"],
            biolink_version,
            kraken_version,
        )

    return integrated_nodes_path, integrated_edges_path


def harmonize_sources(
    sources_config: dict,
    biolink_version: str,
    base_path: Path,
    harmonized_dir_path: Path,
    create_metagraphs: bool,
    metagraph_dir_path: Path,
    zip_inputs_after: bool,
):
    """Harmonize each source that needs it"""
    biolink_client = BiolinkClient(biolink_version)
    for source_name, source_config in sources_config.items():

        harmonize_source(
            source_name,
            source_config,
            biolink_client,
            base_path,
            create_metagraphs,
            harmonized_dir_path,
            metagraph_dir_path,
            zip_inputs_after,
        )


def harmonize_source(
    source_name: str,
    source_config: dict,
    biolink_client: BiolinkClient,
    base_path: Path,
    create_metagraph: bool,
    harmonized_dir_path: Path,
    metagraph_dir_path: Path,
    zip_inputs_after: bool,
):
    """Harmonize a single source to Biolink schema"""
    logging.info(f"Harmonizing {source_name}...")

    # Get output paths
    nodes_output, edges_output = get_harmonized_file_paths(source_name, harmonized_dir_path)

    # Create output directory if it doesn't exist
    nodes_output.parent.mkdir(parents=True, exist_ok=True)
    edges_output.parent.mkdir(parents=True, exist_ok=True)

    # Run source-specific harmonizer
    harmonizers = {
        "kg2": KG2Harmonizer,
        "robokop": RobokopHarmonizer,
        "molepro": MoleProHarmonizer,
        "microbiome-kg": MicrobiomeKGHarmonizer,
        "multiomics-kg": MultiomicsKGHarmonizer,
        "spoke": SpokeHarmonizer,
        "umls": UMLSHarmonizer,
        "lipidmaps": LipidMapsHarmonizer,
        "refmet": RefMetHarmonizer,
        "clingen": ClinGenHarmonizer,
    }

    # Instantiate our harmonizer
    harmonizer = harmonizers[source_name](biolink_client)

    # Construct full input file paths using base path as applicable
    full_input_file_paths = {
        file_slot: base_path / Path(source_config[file_slot])
        for file_slot in ["input_file", "nodes_input", "edges_input"]
        if source_config.get(file_slot)
    }
    logging.info(f"Full input file paths are: {full_input_file_paths}")

    # Unzip input files as needed
    unzip_files(list(full_input_file_paths.values()))

    # Harmonize input files
    if full_input_file_paths.get("input_file"):
        # Note: These are not compatible with the BaseHarmonizer - have separate harmonize() implementations
        harmonizer.harmonize(full_input_file_paths["input_file"], nodes_output, edges_output)
    elif full_input_file_paths.get("nodes_input") and full_input_file_paths.get("edges_input"):
        # Note: These (with split nodes/edges files) use the BaseHarmonizer
        harmonizer.harmonize(
            full_input_file_paths["nodes_input"], full_input_file_paths["edges_input"], nodes_output, edges_output
        )
    else:
        raise ValueError(f"Unknown source type: {source_name}")

    if zip_inputs_after:
        # Zip input files back up
        zip_files(list(full_input_file_paths.values()))

    if create_metagraph:
        generate_metagraph_for_source(nodes_output, edges_output, metagraph_dir_path / source_name, source_name)


def post_process_integrated_kg(
    integrated_nodes_path: Path,
    integrated_edges_path: Path,
    base_path: Path,
    config: dict,
    biolink_version: str,
    kraken_version: str,
):
    """Run all post-processing steps on the unified KG"""
    logging.info("Starting post-processing...")

    if config.get("test_export"):
        logging.info("Generating test files for this kraken build..")
        test_export_config = config["test_export"]
        output_dir = base_path / Path(test_export_config["output_directory"])
        create_test_kg_files(
            integrated_nodes_path, integrated_edges_path, output_dir, num_edges=test_export_config["num_edges"]
        )

    logging.info("Post-processing complete")
