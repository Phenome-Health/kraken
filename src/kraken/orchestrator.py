"""
Main orchestration functions for KRAKEN build
"""

import importlib.metadata
import logging
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

from kraken.biolink_client import BiolinkClient
from kraken.config import KrakenConfig
from kraken.harmonizers.cdes import CDEHarmonizer
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
from kraken.integrate import integrate_sources
from kraken.metagraph import generate_metagraph_for_source
from kraken.post_processing.test_file_generator import create_test_kg_files
from kraken.utils.constants import PROJECT_ROOT
from kraken.utils.kg_io import form_tarball, unzip_files, zip_files
from kraken.utils.logging_config import setup_logging
from kraken.validator import KrakenValidator


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
        "cdes": CDEHarmonizer,
    }

    def __init__(self):
        config_path = Path(f"{PROJECT_ROOT}/config/build_config.yaml")
        with open(config_path) as f:
            config_dict = yaml.safe_load(f)

        self.config = KrakenConfig(**config_dict)
        setup_logging(self.config.log_level)

        self.biolink_client = BiolinkClient(self.config.biolink_version)
        self.validator = KrakenValidator(self.biolink_client)

    def run(self) -> tuple[Path, Path]:
        """Main entry point for building the KRAKEN"""
        logging.info("Starting KRAKEN build...")
        start = time.time()
        logging.info(f"Will include {len(self.config.sources_to_use)} sources: {self.config.sources_to_use}")

        if self.config.steps.harmonize:
            self._harmonize_sources()

        if self.config.steps.integrate:
            self._integrate_sources()

        if self.config.steps.postprocess:
            self._post_process()

        elapsed = time.time() - start
        logging.info(f"Build complete! Took {round(elapsed / 60)} minutes.")

        self._write_build_info(elapsed)

        return self.config.integrated_nodes_path, self.config.integrated_edges_path

    def _write_build_info(self, elapsed_seconds: float) -> None:
        """Write build_info.json to integrated_dir for downstream consumers (e.g. Kestrel /health)."""
        from kraken.build_info import BuildInfo, StepsRun

        try:
            git_commit = (
                subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, cwd=PROJECT_ROOT)
                .decode()
                .strip()
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            git_commit = "unknown"

        try:
            kraken_package_version = importlib.metadata.version("kraken")
        except importlib.metadata.PackageNotFoundError:
            kraken_package_version = "unknown"

        build_info = BuildInfo(
            kg_version=self.config.kraken_version,
            kraken_package_version=kraken_package_version,
            biolink_version=self.config.biolink_version,
            build_timestamp=datetime.now(timezone.utc).isoformat(),
            git_commit=git_commit,
            sources=sorted(self.config.sources_to_use),  # set -> sorted list (JSON-safe, deterministic)
            steps_run=StepsRun(
                harmonize=self.config.steps.harmonize,
                integrate=self.config.steps.integrate,
                postprocess=self.config.steps.postprocess,
            ),
            build_duration_minutes=round(elapsed_seconds / 60, 1),
            kg_label=getattr(self.config, "kg_label", None),
        )

        self.config.integrated_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.integrated_dir / "build_info.json"
        output_path.write_text(build_info.model_dump_json(indent=2))
        logging.info(f"Build info written to {output_path}")

    def _harmonize_sources(self):
        """Harmonize all sources to KRAKEN's Biolink-style semantic layer/schema"""
        logging.info("-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        for source_name in self.config.sources_to_use:
            self._harmonize_source(source_name)

    def _harmonize_source(self, source_name: str):
        """Harmonize a single source to Biolink schema"""
        logging.info(f"Harmonizing {source_name}...")
        if self.config.options.validation_only:
            logging.warning("Skipping harmonization step because you selected validation_only=True in build config")

        source_config = self.config.sources[source_name]

        # Get paths
        nodes_output, edges_output = self.config.all_harmonized_paths_resolved[source_name]

        # Create output directory if it doesn't exist
        nodes_output.parent.mkdir(parents=True, exist_ok=True)

        # Instantiate our harmonizer
        harmonizer = self.HARMONIZERS[source_name](self.biolink_client)

        if not self.config.options.validation_only:
            # Unzip input files as needed
            unzip_files(self.config.all_source_input_paths_resolved[source_name])

            harmonizer.harmonize(
                nodes_output=nodes_output,
                edges_output=edges_output,
                input_file=source_config.input_file_resolved,
                nodes_input=source_config.nodes_input_resolved,
                edges_input=source_config.edges_input_resolved,
            )

            if self.config.zip_inputs_after:
                zip_files(self.config.all_source_input_paths_resolved[source_name])

        if self.config.options.validate_output or self.config.options.validation_only:
            self.validator.validate(nodes_output, edges_output, harmonizer.source_infores)

        if self.config.create_metagraphs and not self.config.options.validation_only:
            generate_metagraph_for_source(
                nodes_path=nodes_output,
                edges_path=edges_output,
                output_dir=self.config.metagraph_dir / source_name,
                graph_name=source_name,
            )

    def _integrate_sources(self):
        """Integrate sources into unified KG with entity resolution"""
        logging.info("-------------------------- INTEGRATING SOURCES -----------------------------------------------")
        if self.config.options.validation_only:
            logging.warning("Skipping integration step because you selected validation_only=True in build config")

        if not self.config.options.validation_only:
            integrate_sources(self.config, self.biolink_client)

        if self.config.options.validate_output or self.config.options.validation_only:
            self.validator.validate(
                self.config.integrated_nodes_path, self.config.integrated_edges_path, integrated=True
            )

        if not self.config.options.validation_only:

            tarball_component_paths = [self.config.integrated_nodes_path, self.config.integrated_edges_path]

            if self.config.create_metagraphs:
                metagraph_path = generate_metagraph_for_source(
                    nodes_path=self.config.integrated_nodes_path,
                    edges_path=self.config.integrated_edges_path,
                    output_dir=self.config.metagraph_dir,
                    graph_name="kraken",
                    graph_version=self.config.kraken_version,
                )
                tarball_component_paths.append(metagraph_path)

            form_tarball(
                tarball_component_paths,
                self.config.integrated_dir,
                tarball_name=f"kraken_{self.config.kraken_version}.tar.gz",
            )

    def _post_process(self):
        """Run all post-processing steps on the unified KG"""
        logging.info("------------------------------ POST-PROCESSING -----------------------------------------------")

        if self.config.post_processing:

            if self.config.post_processing.test_export:
                logging.info("Generating test files for this kraken build..")
                test_export_config = self.config.post_processing.test_export

                test_nodes_path, test_edges_path = create_test_kg_files(
                    nodes_path=self.config.integrated_nodes_path,
                    edges_path=self.config.integrated_edges_path,
                    output_dir=self.config.test_export_dir,
                    num_edges=test_export_config.num_edges,
                )
                test_tarball_component_paths = [test_nodes_path, test_edges_path]

                if self.config.create_metagraphs:
                    metagraph_path = generate_metagraph_for_source(
                        nodes_path=test_nodes_path,
                        edges_path=test_edges_path,
                        output_dir=self.config.metagraph_dir,
                        graph_name="kraken_test",
                    )
                    test_tarball_component_paths.append(metagraph_path)

                form_tarball(
                    test_tarball_component_paths,
                    self.config.integrated_dir,
                    tarball_name=f"kraken_{self.config.kraken_version}_test.tar.gz",
                )

        logging.info("Post-processing complete!")
