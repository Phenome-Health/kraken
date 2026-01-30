"""Validator for KRAKEN-ified JSONL files."""

import logging
import re
from pathlib import Path

from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.constants import (
    AGENT_TYPE,
    CATEGORIES,
    CORE_NODE_PROPERTIES,
    EDGE_PROPERTIES,
    EQUIVALENT_IDS,
    ID,
    KNOWLEDGE_LEVEL,
    OBJECT,
    PREDICATE,
    PRIMARY_KS,
    PROVIDED_BY,
    SUBJECT,
    TYPE,
)
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl


class KrakenValidator:
    """Validates KRAKEN-harmonized JSONL node and edge files."""

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client

    def validate(self, nodes_path: Path, edges_path: Path, source_infores: str | None = None) -> None:
        """
        Validate KRAKEN-harmonized node and edge files.

        Raises:
            ValueError: If any validation check fails.
        """
        logging.info(f"Validating {nodes_path} and {edges_path}...")

        self.validate_nodes(nodes_path, source_infores)
        self.validate_edges(edges_path)

        # TODO: attributes (nodes and edges), types within lists, non-core edge props, aggregator ks.. qualifiers?
        # TODO: name, synonyms, none strings?

        logging.info("Validation complete! All checks passed.")

    def validate_nodes(self, nodes_path: Path, source_infores: str | None = None) -> None:
        logging.info(f"Starting to validate {nodes_path}")

        self._check_file_suffix(nodes_path)

        required_fields = [ID, CATEGORIES, PROVIDED_BY, EQUIVALENT_IDS]

        for node in stream_nodes_from_jsonl(nodes_path):

            # Required fields present and non-empty
            for field in required_fields:
                if field not in node:
                    self._raise(f"Missing required field node.{field}", node)
                if not node[field]:
                    self._raise(f"Field node.{field} cannot be empty", node)

            # Type checks
            for prop, expected_type in CORE_NODE_PROPERTIES.items():
                if prop in node:
                    value = node[prop]
                    if not isinstance(value, expected_type):
                        self._raise(
                            f"Field node.{prop} has type {type(value).__name__}, expected {expected_type.__name__}",
                            node,
                        )

            # ID must appear in equivalent_ids
            if node[ID] not in node[EQUIVALENT_IDS]:
                self._raise(f"node.{ID} must appear in node.{EQUIVALENT_IDS}", node)

            # Source provenance check (if source_infores provided)
            if source_infores and source_infores not in node[PROVIDED_BY]:
                self._raise(f"node.{PROVIDED_BY} must contain '{source_infores}'", node)

            # Categories must be valid Biolink categories
            for category in node[CATEGORIES]:
                if not self.is_valid_category_format(category):
                    self._raise(f"Category '{category}' is not in proper biolink:PascalCase format", node)
                if not (category in self.biolink.categories or self.biolink.toolkit.is_mixin(category)):
                    # NOTE: biolink has bug where some category mixins are not descendants of NamedThing
                    self._raise(f"Category '{category}' does not exist in Biolink", node)

    def validate_edges(self, edges_path: Path) -> None:
        logging.info(f"Starting to validate {edges_path}")

        self._check_file_suffix(edges_path)

        required_fields = [SUBJECT, PREDICATE, OBJECT, PRIMARY_KS, KNOWLEDGE_LEVEL, AGENT_TYPE]

        for edge in stream_edges_from_jsonl(edges_path):
            # Required fields present and non-empty
            for field in required_fields:
                if field not in edge:
                    self._raise(f"Edge missing required field '{field}'", edge)
                if not edge[field]:
                    self._raise(f"Field edge.{field} cannot be empty", edge)

            # Type checks
            for prop, info in EDGE_PROPERTIES.items():
                if prop in edge:
                    value = edge[prop]
                    expected_type = info[TYPE]
                    if not isinstance(value, expected_type):
                        self._raise(
                            f"Field edge.{prop} has type {type(value).__name__}, expected {expected_type.__name__}",
                            edge,
                        )

            # Predicate must be valid Biolink predicate
            predicate = edge[PREDICATE]
            if not self.is_valid_predicate_format(predicate):
                self._raise(f"Predicate '{predicate}' is not in proper biolink:snake_case format", edge)
            if predicate not in self.biolink.predicates:
                self._raise(f"Invalid predicate '{predicate}'", edge)

            # Knowledge level must be valid per Biolink
            if edge[KNOWLEDGE_LEVEL] not in self.biolink.knowledge_levels:
                self._raise(
                    f"Invalid edge.knowledge_level '{edge[KNOWLEDGE_LEVEL]}'. "
                    f"Valid options are: {self.biolink.knowledge_levels}",
                    edge,
                )
            # Agent type must be valid per Biolink
            if edge[AGENT_TYPE] not in self.biolink.agent_types:
                self._raise(
                    f"Invalid edge.agent_type '{edge[AGENT_TYPE]}'. Valid options are: {self.biolink.agent_types}",
                    edge,
                )

        logging.info(f"Edge validation complete for {edges_path}")

    def _check_file_suffix(self, path: Path) -> None:
        if path.suffix != ".jsonl":
            self._raise(f"File must have .jsonl suffix: {path}")

    @staticmethod
    def _raise(message: str, item: dict | None = None):
        item_str = f"Full item: {item}" if item else None
        raise ValueError(f"{message}. {item_str}")

    @staticmethod
    def is_valid_category_format(category: str) -> bool:
        """Check that category is in 'biolink:PascalCase' format."""
        return bool(re.match(r"^biolink:[A-Z][a-zA-Z]*$", category))

    @staticmethod
    def is_valid_predicate_format(predicate: str) -> bool:
        """Check that predicate is in 'biolink:snake_case' format."""
        return bool(re.match(r"^biolink:[a-z][a-z_]*$", predicate))
