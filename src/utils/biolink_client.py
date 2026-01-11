# biolink_client.py
from typing import Set, List
import logging

from .general import to_list

from bmt import Toolkit


class BiolinkClient:
    """Client for Biolink Model operations"""

    def __init__(self, biolink_version: str):
        biolink_url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{biolink_version}/biolink-model.yaml"
        logging.info(f"Initializing Biolink Model Toolkit for version {biolink_version}...")
        self.toolkit = Toolkit(schema=biolink_url)
        self.version = biolink_version

    def filter_to_leaf_categories(self, categories: str | list[str] | set[str]) -> list[str]:
        """Remove ancestral categories, keeping only the most specific (leaf) ones"""
        categories = set(to_list(categories))
        all_proper_ancestors = set()

        for category in categories:
            proper_ancestors = set(self.toolkit.get_ancestors(
                category,
                formatted=True,
                mixin=True,
                reflexive=False
            ))
            all_proper_ancestors |= proper_ancestors

        return list(categories - all_proper_ancestors)
