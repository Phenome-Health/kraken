#!/usr/bin/env python3
"""
KRAKEN Build System
Main entry point for building the unified knowledge graph
"""

import logging
from pathlib import Path
import time
import yaml

from src.orchestrator import run_kg_build
from src.utils.logging_config import setup_logging


def main():
    setup_logging()

    config_path = Path("config/build_config.yaml")
    with open(config_path) as f:
        config = yaml.safe_load(f)

    logging.info("Starting KRAKEN build...")
    start = time.time()
    final_kg_path = run_kg_build(config)
    logging.info(f"Build complete! Took {round((time.time() - start) / 60)} minutes. KRAKEN saved to: {final_kg_path}")


if __name__ == "__main__":
    main()