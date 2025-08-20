"""
Main orchestration functions for PhenomeKG build
Functional approach to workflow management
"""

from pathlib import Path
from typing import Dict, List
import logging

from .harmonizers import get_harmonizer
from .integration.entity_resolution import integrate_sources
from .post_processing.arango_export import prepare_for_arango
from .post_processing.biomapper_export import export_for_biomapper
from .utils.kg_io import load_kg, save_kg


def run_kg_build(config: dict) -> Path:
    """Main orchestration function for building PhenomeKG"""
    logging.info("Starting PhenomeKG build...")

    # Phase 1: Harmonize all sources to Biolink schema
    harmonized_artifacts = harmonize_all_sources(config['sources'])

    # Phase 2: Integrate into unified KG with entity resolution
    phenomekg_path = integrate_sources(harmonized_artifacts, config['integration'])

    # Phase 3: Post-processing steps
    if 'post_processing' in config:
        post_process_unified_kg(phenomekg_path, config['post_processing'])

    logging.info(f"Build complete: {phenomekg_path}")
    return phenomekg_path


def harmonize_all_sources(sources_config: dict) -> Dict[str, Path]:
    """Harmonize each source that needs it"""
    harmonized_paths = {}

    for source_name, source_config in sources_config.items():
        if needs_harmonization(source_name, source_config):
            harmonized_paths[source_name] = harmonize_source(source_name, source_config)
        else:
            logging.info(f"Skipping {source_name} - already harmonized")
            harmonized_paths[source_name] = Path(source_config['harmonized_output'])

    return harmonized_paths


def harmonize_source(source_name: str, config: dict) -> Path:
    """Harmonize a single source to Biolink schema"""
    input_path = Path(config['input_path'])
    output_path = Path(config['harmonized_output'])

    logging.info(f"Harmonizing {source_name}...")

    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get source-specific harmonizer and run it
    harmonizer = get_harmonizer(source_name)
    harmonizer(input_path, output_path, config.get('harmonization_rules', {}))

    logging.info(f"Harmonized {source_name} -> {output_path}")
    return output_path


def needs_harmonization(source_name: str, config: dict) -> bool:
    """Check if harmonization step needs to run based on file timestamps"""
    input_path = Path(config['input_path'])
    output_path = Path(config['harmonized_output'])

    if not output_path.exists():
        logging.info(f"{source_name} needs harmonization: output doesn't exist")
        return True

    if not input_path.exists():
        logging.warning(f"{source_name} input file doesn't exist: {input_path}")
        return False

    input_mtime = input_path.stat().st_mtime
    output_mtime = output_path.stat().st_mtime

    needs_update = input_mtime > output_mtime
    if needs_update:
        logging.info(f"{source_name} needs harmonization: input is newer")

    return needs_update


def post_process_unified_kg(unified_kg_path: Path, config: dict):
    """Run all post-processing steps on the unified KG"""
    logging.info("Starting post-processing...")

    # Load the unified KG once
    unified_kg = load_kg(unified_kg_path)

    # Step 1: Prepare ArangoDB version
    if 'arango_export' in config:
        arango_config = config['arango_export']
        if needs_arango_export(unified_kg_path, arango_config):
            logging.info("Preparing ArangoDB export...")
            prepare_for_arango(unified_kg, arango_config)
        else:
            logging.info("Skipping ArangoDB export - up to date")

    # Step 2: Export for biomapper
    if 'biomapper_export' in config:
        biomapper_config = config['biomapper_export']
        if needs_biomapper_export(unified_kg_path, biomapper_config):
            logging.info("Exporting for biomapper...")
            export_for_biomapper(unified_kg, biomapper_config)
        else:
            logging.info("Skipping biomapper export - up to date")

    logging.info("Post-processing complete")


def needs_arango_export(unified_kg_path: Path, config: dict) -> bool:
    """Check if ArangoDB export needs to run"""
    output_path = Path(config['output_path'])

    if not output_path.exists():
        return True

    return unified_kg_path.stat().st_mtime > output_path.stat().st_mtime


def needs_biomapper_export(unified_kg_path: Path, config: dict) -> bool:
    """Check if biomapper export needs to run"""
    output_dir = Path(config['output_directory'])

    if not output_dir.exists():
        return True

    # Check if any of the expected output files are missing or outdated
    unified_mtime = unified_kg_path.stat().st_mtime

    for file_path in output_dir.iterdir():
        if file_path.stat().st_mtime < unified_mtime:
            return True

    return False