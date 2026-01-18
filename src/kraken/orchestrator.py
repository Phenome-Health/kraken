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
from kraken.utils.kg_io import form_tarball, get_harmonized_file_paths, unzip_files, zip_files
from kraken.utils.metagraph import generate_metagraph_for_source


class KrakenBuildOrchestrator:
    """Main orchestrator for building the KRAKEN knowledge graph"""

    HARMONIZERS = {
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

    def __init__(self, config: dict):
        self.config = config
        self.biolink_version = config["biolink_version"]
        self.kraken_version = config["kraken_version"]
        self.base_path = Path(config["base_path"]) if config.get("base_path") else PROJECT_ROOT
        self.biolink_client: BiolinkClient = BiolinkClient(self.biolink_version)

        # Paths
        self.harmonized_dir = self.base_path / config["harmonization"]["output_directory"]
        self.metagraph_dir = self.base_path / config["metagraph"]["output_directory"]
        self.integrated_dir = self.base_path / config["integration"]["output_directory"]
        self.integrated_nodes_path = self.integrated_dir / f"kraken_nodes_{self.kraken_version}.jsonl"
        self.integrated_edges_path = self.integrated_dir / f"kraken_edges_{self.kraken_version}.jsonl"

        # Options
        self.create_metagraphs = bool(config["options"].get("metagraph_creation"))
        self.zip_inputs_after = bool(config["harmonization"].get("zip_inputs_after"))
        self.sources_to_use: set[str] = self._resolve_sources()

    def run(self) -> tuple[Path, Path]:
        """Main entry point for building the KRAKEN"""
        if self.config["steps"].get("harmonize"):
            self._harmonize_sources()

        if self.config["steps"].get("integrate"):
            self._integrate_sources()

        if self.config["steps"].get("postprocess"):
            self._post_process()

        return self.integrated_nodes_path, self.integrated_edges_path

    def _resolve_sources(self) -> set[str]:
        """Figure out which sources to use based on the build config settings"""
        all_sources = set(self.config["sources"])
        include_sources = set(to_list(self.config["options"].get("include_sources")))
        exclude_sources = set(to_list(self.config["options"].get("exclude_sources")))

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
        return sources_to_use

    def _harmonize_sources(self):
        """Harmonize all sources to KRAKEN's Biolink-style semantic layer/schema"""
        logging.info("-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        for source_name in self.sources_to_use:
            self._harmonize_source(source_name)

    def _harmonize_source(self, source_name: str):
        """Harmonize a single source to Biolink schema"""
        logging.info(f"Harmonizing {source_name}...")

        # Get output paths
        nodes_output, edges_output = get_harmonized_file_paths(source_name, self.harmonized_dir)

        # Create output directory if it doesn't exist
        nodes_output.parent.mkdir(parents=True, exist_ok=True)
        edges_output.parent.mkdir(parents=True, exist_ok=True)

        # Instantiate our harmonizer
        harmonizer = self.HARMONIZERS[source_name](self.biolink_client)

        # Construct full input file paths using base path as applicable
        source_config = self.config["sources"][source_name]
        full_input_file_paths = {
            file_slot: self.base_path / Path(source_config[file_slot])
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

        if self.zip_inputs_after:
            zip_files(list(full_input_file_paths.values()))

        if self.create_metagraphs:
            generate_metagraph_for_source(
                nodes_path=nodes_output,
                edges_path=edges_output,
                output_dir=self.metagraph_dir / source_name,
                graph_name=source_name,
            )

    def _integrate_sources(self):
        """Integrate sources into unified KG with entity resolution"""
        logging.info("-------------------------- INTEGRATING SOURCES -----------------------------------------------")
        integrate_sources(
            source_names=self.sources_to_use,
            integrated_nodes_path=self.integrated_nodes_path,
            integrated_edges_path=self.integrated_edges_path,
            harmonized_dir_path=self.harmonized_dir,
            config=self.config,
        )
        tarball_component_paths = [self.integrated_nodes_path, self.integrated_edges_path]

        if self.create_metagraphs:
            metagraph_path = generate_metagraph_for_source(
                nodes_path=self.integrated_nodes_path,
                edges_path=self.integrated_edges_path,
                output_dir=self.metagraph_dir,
                graph_name="kraken",
            )
            tarball_component_paths.append(metagraph_path)

        form_tarball(tarball_component_paths, self.integrated_dir)

    def _post_process(self):
        """Run all post-processing steps on the unified KG"""
        logging.info("------------------------------ POST-PROCESSING -----------------------------------------------")
        post_config = self.config["post_processing"]

        if post_config.get("test_export"):
            logging.info("Generating test files for this kraken build..")
            test_export_config = post_config["test_export"]
            test_nodes_path, test_edges_path = create_test_kg_files(
                nodes_path=self.integrated_nodes_path,
                edges_path=self.integrated_edges_path,
                output_dir=self.base_path / Path(test_export_config["output_directory"]),
                num_edges=test_export_config["num_edges"],
            )
            tarball_component_paths = [test_nodes_path, test_edges_path]

            if self.create_metagraphs:
                metagraph_path = generate_metagraph_for_source(
                    nodes_path=test_nodes_path,
                    edges_path=test_edges_path,
                    output_dir=self.metagraph_dir,
                    graph_name="kraken_test",
                )
                tarball_component_paths.append(metagraph_path)

            form_tarball(tarball_component_paths, self.integrated_dir)

        logging.info("Post-processing complete")


def run_build(config: dict) -> tuple[Path, Path]:
    """Main orchestration function for building the KRAKEN"""
    return KrakenBuildOrchestrator(config).run()
