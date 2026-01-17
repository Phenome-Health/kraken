# biolink_client.py
import logging

from kraken.utils.general import to_list

from bmt import Toolkit


class BiolinkClient:
    """Client for Biolink Model operations"""

    def __init__(self, biolink_version: str):
        biolink_url = f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{biolink_version}/biolink-model.yaml"
        logging.info(f"Initializing Biolink Model Toolkit for version {biolink_version}...")
        self.toolkit = Toolkit(schema=biolink_url)
        self.version = biolink_version
        self.all_categories = set(self.toolkit.get_descendants("biolink:NamedThing",
                                                               formatted=True,
                                                               mixin=True,
                                                               reflexive=True))
        self.all_predicates = set(self.toolkit.get_descendants("biolink:related_to",
                                                               formatted=True,
                                                               mixin=True,
                                                               reflexive=True))

    def filter_to_leaf_categories(self, categories: str | list[str] | set[str]) -> list[str]:
        """Remove ancestral categories, keeping only the most specific (leaf) categories"""
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


def main():
    bc = BiolinkClient("4.2.5")
    print(bc.all_predicates)
    print(bc.all_predicates)
    print(bc.toolkit.get_permissible_value_ids_for_slot("knowledge_level"))