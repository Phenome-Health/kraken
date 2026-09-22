"""
Main orchestration functions for KRAKEN build
"""

import importlib.metadata
import logging
import subprocess
import sys
import time
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path

import yaml

from kraken.biolink_client import BiolinkClient
from kraken.config import KrakenConfig
from kraken.harmonizers.babel import BabelHarmonizer
from kraken.harmonizers.bio_age import BioAgeHarmonizer
from kraken.harmonizers.bio_bmi import BioBMIHarmonizer
from kraken.harmonizers.cdes import CDEHarmonizer
from kraken.harmonizers.clingen import ClinGenHarmonizer
from kraken.harmonizers.ctkg import CTKGHarmonizer
from kraken.harmonizers.dakg import DAKGHarmonizer
from kraken.harmonizers.kg2 import KG2Harmonizer
from kraken.harmonizers.lipidmaps import LipidMapsHarmonizer
from kraken.harmonizers.loinc import LoincHarmonizer
from kraken.harmonizers.long_covid import LongCovidHarmonizer
from kraken.harmonizers.microbiome_kg import MicrobiomeKGHarmonizer
from kraken.harmonizers.multiomics_kg import MultiomicsKGHarmonizer
from kraken.harmonizers.ncbigene import NCBIGeneHarmonizer
from kraken.harmonizers.pgs_catalog import PGSCatalogHarmonizer
from kraken.harmonizers.refmet import RefMetHarmonizer
from kraken.harmonizers.robokop import RobokopHarmonizer
from kraken.harmonizers.translator_kg_open import TranslatorKGOpenHarmonizer
from kraken.harmonizers.umls import UMLSHarmonizer
from kraken.integrate import integrate_sources
from kraken.metagraph import generate_metagraph_for_source
from kraken.post_processing.test_file_generator import create_test_kg_files
from kraken.utils.constants import PROJECT_ROOT
from kraken.utils.kg_io import unzip_files, zip_files
from kraken.utils.logging_config import setup_logging
from kraken.validator import KrakenValidator


def harmonization_order(sources: Iterable[str]) -> list[str]:
    """The order to harmonize sources in: alphabetical, except Babel last, since it reads the others' output."""
    return sorted(sources, key=lambda source: (source == "babel", source))


class KrakenBuildOrchestrator:
    """Main orchestrator for building the KRAKEN knowledge graph"""

    HARMONIZERS = {
        "kg2": KG2Harmonizer,
        "robokop": RobokopHarmonizer,
        "ctkg": CTKGHarmonizer,
        "dakg": DAKGHarmonizer,
        "microbiome-kg": MicrobiomeKGHarmonizer,
        "ncbigene": NCBIGeneHarmonizer,
        "multiomics-kg": MultiomicsKGHarmonizer,
        "umls": UMLSHarmonizer,
        "lipidmaps": LipidMapsHarmonizer,
        "loinc": LoincHarmonizer,
        "refmet": RefMetHarmonizer,
        "clingen": ClinGenHarmonizer,
        "cdes": CDEHarmonizer,
        "translator-kg-open": TranslatorKGOpenHarmonizer,
        "pgs-catalog": PGSCatalogHarmonizer,
        "bio-bmi": BioBMIHarmonizer,
        "bio-age": BioAgeHarmonizer,
        "long-covid": LongCovidHarmonizer,
        "babel": BabelHarmonizer,
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
        self._log_source_versions()
        self._confirm_source_versions()

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

    def _log_source_versions(self) -> None:
        """Log each in-use source's declared version next to its resolved input path(s), so any drift
        between the configured version and the actual input files is easy to eyeball before a build."""
        logging.info("Source versions for this build (version <- input path):")
        for source_name in sorted(self.config.sources_to_use):
            version = self.config.sources[source_name].version
            paths = ", ".join(str(p) for p in self.config.all_source_input_paths_resolved[source_name])
            logging.info(f"  {source_name}: {version!r} <- {paths}")

    def _confirm_source_versions(self) -> None:
        """Ask the user to confirm the listed source versions before the build proceeds. Skipped when
        disabled (options.confirm_source_versions=False) or when stdin is not interactive (cron/CI), so
        it never blocks non-interactive runs."""
        if not self.config.options.confirm_source_versions:
            return
        if not sys.stdin.isatty():
            logging.info("stdin is not interactive; skipping source-version confirmation.")
            return
        try:
            response = input("Proceed with the build using these source versions? [y/N]: ").strip().lower()
        except EOFError:
            response = ""
        if response not in {"y", "yes"}:
            logging.info("Build aborted: source versions not confirmed.")
            raise SystemExit(0)

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
            source_versions=getattr(self.config, "source_versions", None),
        )

        self.config.integrated_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.integrated_dir / "build_info.json"
        output_path.write_text(build_info.model_dump_json(indent=2))
        logging.info(f"Build info written to {output_path}")

    def _harmonize_sources(self):
        """Harmonize all sources to KRAKEN's Biolink-style semantic layer/schema"""
        logging.info("-------------------------- HARMONIZING SOURCES -----------------------------------------------")
        for source_name in harmonization_order(self.config.sources_to_use):
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

        # Instantiate our harmonizer (its provenance id comes from build_config: sources.<name>.source_id),
        # telling it which other sources' edges to drop because we ingest those directly
        # (build_config: sources.<other>.drop_from_other_sources).
        extra_arguments = {}
        if source_name == "babel":
            # Babel decides which structure-only cliques to keep by what the other sources reference, so it reads
            # their harmonized nodes (and is harmonized last -- see harmonization_order).
            extra_arguments["other_sources_nodes"] = self.config.harmonized_nodes_paths_of_build_sources(
                other_than=source_name
            )
        harmonizer = self.HARMONIZERS[source_name](
            self.biolink_client,
            source_id=source_config.source_id,
            auto_source_exclusions=self.config.auto_source_exclusions(source_name),
            **extra_arguments,
        )

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
            # Report the curies biomapper2 couldn't fully normalize. Driven from here, not from the
            # harmonizer, because the single-file harmonizers override harmonize() -- so this is the
            # one place that runs for every source.
            harmonizer.log_normalization_report()
            harmonizer.log_taxon_report()

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
                source_versions={source_name: self.config.sources[source_name].version},
                biolink_version=self.config.biolink_version,
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
            if self.config.create_metagraphs:
                generate_metagraph_for_source(
                    nodes_path=self.config.integrated_nodes_path,
                    edges_path=self.config.integrated_edges_path,
                    output_dir=self.config.metagraph_dir,
                    graph_name="kraken",
                    graph_version=self.config.kraken_version,
                    source_versions=self.config.source_versions,
                    biolink_version=self.config.biolink_version,
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

                if self.config.create_metagraphs:
                    generate_metagraph_for_source(
                        nodes_path=test_nodes_path,
                        edges_path=test_edges_path,
                        output_dir=self.config.metagraph_dir,
                        graph_name="kraken_test",
                        graph_version=self.config.kraken_version,
                        source_versions=self.config.source_versions,
                        biolink_version=self.config.biolink_version,
                    )

        logging.info("Post-processing complete!")
