"""
Logging configuration for PhenomeKG build
"""

import logging
import sys


def setup_logging(level=logging.INFO):
    """Setup logging configuration"""

    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    # Setup console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)

    # Setup file handler
    file_handler = logging.FileHandler('kraken_build.log')
    file_handler.setFormatter(formatter)

    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.addHandler(console_handler)
    root_logger.addHandler(file_handler)

    # Prevent duplicate logs
    root_logger.propagate = False
