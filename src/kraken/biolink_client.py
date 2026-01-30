# biolink_client.py
import logging

from bmt import Toolkit

from kraken.utils.constants import ROOT_CATEGORY, ROOT_PREDICATE
from kraken.utils.general import to_list


class BiolinkClient:
    """Client for Biolink Model operations"""

    def __init__(self, biolink_version: str):
        biolink_url = (
            f"https://raw.githubusercontent.com/biolink/biolink-model/refs/tags/v{biolink_version}/biolink-model.yaml"
        )
        logging.info(f"Initializing Biolink Model Toolkit for version {biolink_version}...")
        self.toolkit = Toolkit(schema=biolink_url)
        self.version = biolink_version
        self.categories = set(self.toolkit.get_descendants(ROOT_CATEGORY, formatted=True, mixin=True, reflexive=True))
        self.predicates = set(self.toolkit.get_descendants(ROOT_PREDICATE, formatted=True, mixin=True, reflexive=True))
        kl_enum = self.toolkit.view.schema.enums.get("KnowledgeLevelEnum")
        self.knowledge_levels = set(kl_enum.permissible_values.keys())
        at_enum = self.toolkit.view.schema.enums.get("AgentTypeEnum")
        self.agent_types = set(at_enum.permissible_values.keys())

    def filter_to_leaf_categories(self, categories: str | list[str] | set[str]) -> list[str]:
        """Remove ancestral categories, keeping only the most specific (leaf) categories"""
        categories = set(to_list(categories))
        all_proper_ancestors = set()

        for category in categories:
            proper_ancestors = set(self.toolkit.get_ancestors(category, formatted=True, mixin=True, reflexive=False))
            all_proper_ancestors |= proper_ancestors

        return list(categories - all_proper_ancestors)


def main():
    bc = BiolinkClient("4.2.5")

    print(bc.categories, "\n")
    print(bc.predicates, "\n")
    print(bc.agent_types, "\n")
    print(bc.knowledge_levels)


if __name__ == "__main__":
    main()
