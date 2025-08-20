#!/usr/bin/env python3
"""
PhenomeKG Build System
Main entry point for building the unified knowledge graph
"""

import logging
from pathlib import Path
import yaml

from src.orchestrator import run_kg_build
from src.utils.logging_config import setup_logging


def main():
    setup_logging()

    config_path = Path("config/build_config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    logging.info("Starting PhenomeKG build...")
    final_kg_path = run_kg_build(config)
    logging.info(f"Build complete! PhenomeKG saved to: {final_kg_path}")


if __name__ == "__main__":
    main()