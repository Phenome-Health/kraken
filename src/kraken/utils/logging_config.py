"""
Logging configuration for KRAKEN build
"""

import logging
import sys

from kraken.utils.constants import PROJECT_ROOT


def setup_logging(level: str = "INFO"):
    """Setup logging configuration"""

    # Clear any existing handlers
    root_logger = logging.getLogger()
    root_logger.handlers = []

    # Create formatter
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    # Setup console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)

    # Setup file handler
    file_handler = logging.FileHandler(PROJECT_ROOT / "kraken_build.log")
    file_handler.setFormatter(formatter)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper()))
    root_logger.addHandler(console_handler)
    root_logger.addHandler(file_handler)

    # Prevent duplicate logs
    root_logger.propagate = False
