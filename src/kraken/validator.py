"""Validator for KRAKEN-ified JSONL files."""

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from kraken.schema import EdgeModel, NodeModel
from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl


@dataclass
class ValidationError:
    """A single validation error."""

    message: str
    item: dict | None = None


@dataclass
class ValidationSummary:
    """Collects and summarizes validation errors."""

    errors_by_type: dict[str, list[ValidationError]] = field(default_factory=lambda: defaultdict(list))
    max_examples_per_type: int = 3

    def add_error(self, error_type: str, message: str, item: dict | None = None):
        self.errors_by_type[error_type].append(ValidationError(message, item))

    @property
    def total_errors(self) -> int:
        return sum(len(errors) for errors in self.errors_by_type.values())

    @property
    def has_errors(self) -> bool:
        return self.total_errors > 0

    def get_summary(self) -> str:
        if not self.has_errors:
            return "No validation errors found."

        lines = [
            f"Validation failed with {self.total_errors} total error(s) across {len(self.errors_by_type)} type(s):\n"
        ]

        for error_type, errors in sorted(self.errors_by_type.items()):
            lines.append(f"  [{error_type}] - {len(errors)} occurrence(s)")
            for error in errors[: self.max_examples_per_type]:
                lines.append(f"    • {error.message}")
                if error.item:
                    item_id = error.item.get("id") or error.item.get("subject", "unknown")
                    lines.append(f"      Item ID: {item_id}")
            if len(errors) > self.max_examples_per_type:
                lines.append(f"    ... and {len(errors) - self.max_examples_per_type} more")

        return "\n".join(lines)


class KrakenValidator:
    """Validates KRAKEN-harmonized JSONL node and edge files."""

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client
        self.summary = ValidationSummary()

    def validate(self, nodes_path: Path, edges_path: Path, source_infores: str | None = None) -> None:
        """
        Validate KRAKEN-harmonized node and edge files.

        Raises:
            ValueError: If any validation check fails.
        """
        logging.info(f"Validating {nodes_path} and {edges_path}...")

        # Reset summary for fresh validation run
        self.summary = ValidationSummary()

        self.validate_nodes(nodes_path, source_infores)
        self.validate_edges(edges_path)

        # TODO: attributes (nodes and edges), types within lists, non-core edge props, aggregator ks.. qualifiers?
        # TODO: name, synonyms, none strings?

        if self.summary.has_errors:
            raise ValueError(self.summary.get_summary())

        logging.info("Validation complete! All checks passed.")

    def validate_nodes(self, nodes_path: Path, source_infores: str | None = None) -> None:
        logging.info(f"Starting to validate {nodes_path}")

        self._check_file_suffix(nodes_path)

        required_props = [p for p in NodeModel.all_properties().values() if p.required]

        for node in stream_nodes_from_jsonl(nodes_path):

            # Required fields present and non-empty
            for prop in required_props:
                if prop.name not in node:
                    self._add_error("missing_required_field", f"Missing required field node.{prop.name}", node)
                elif not node[prop.name]:
                    self._add_error("empty_required_field", f"Field node.{prop.name} cannot be empty", node)

            # Type checks
            for prop in NodeModel.all_properties().values():
                if prop.name in node:
                    value = node[prop.name]
                    if not isinstance(value, prop.type):
                        self._add_error(
                            "invalid_type",
                            f"Field node.{prop.name} has type {type(value).__name__}, expected {prop.type.__name__}",
                            node,
                        )

            # ID must appear in equivalent_ids
            if NodeModel.equivalent_ids.name in node and NodeModel.id.name in node:
                if node[NodeModel.id.name] not in node[NodeModel.equivalent_ids.name]:
                    self._add_error(
                        "id_not_in_equivalent_ids",
                        f"node.{NodeModel.id.name} must appear in node.{NodeModel.equivalent_ids.name}",
                        node,
                    )

            # Source provenance check (if source_infores provided)
            if source_infores and NodeModel.provided_by.name in node:
                if source_infores not in node[NodeModel.provided_by.name]:
                    self._add_error(
                        "missing_source_provenance",
                        f"node.{NodeModel.provided_by.name} must contain '{source_infores}'",
                        node,
                    )

            # Categories must be valid Biolink categories
            if NodeModel.categories.name in node:
                for category in node[NodeModel.categories.name]:
                    if not self.is_valid_category_format(category):
                        self._add_error(
                            "invalid_category_format",
                            f"Category '{category}' is not in proper biolink:PascalCase format",
                            node,
                        )
                    elif not (category in self.biolink.categories or self.biolink.toolkit.is_mixin(category)):
                        # NOTE: biolink has bug where some category mixins are not descendants of NamedThing
                        self._add_error(
                            "invalid_category",
                            f"Category '{category}' does not exist in Biolink",
                            node,
                        )

    def validate_edges(self, edges_path: Path) -> None:
        logging.info(f"Starting to validate {edges_path}")

        self._check_file_suffix(edges_path)

        required_props = [p for p in EdgeModel.all_properties().values() if p.required]

        for edge in stream_edges_from_jsonl(edges_path):
            # Required fields present and non-empty
            for prop in required_props:
                if prop.name not in edge:
                    self._add_error("missing_required_field", f"Edge missing required field '{prop.name}'", edge)
                elif not edge[prop.name]:
                    self._add_error("empty_required_field", f"Field edge.{prop.name} cannot be empty", edge)

            # Type checks
            for prop in EdgeModel.all_properties().values():
                if prop.name in edge:
                    value = edge[prop.name]
                    if not isinstance(value, prop.type):
                        self._add_error(
                            "invalid_type",
                            f"Field edge.{prop.name} has type {type(value).__name__}, expected {prop.type.__name__}",
                            edge,
                        )

            # Predicate must be valid Biolink predicate
            if EdgeModel.predicate.name in edge:
                predicate = edge[EdgeModel.predicate.name]
                if not self.is_valid_predicate_format(predicate):
                    self._add_error(
                        "invalid_predicate_format",
                        f"Predicate '{predicate}' is not in proper biolink:snake_case format",
                        edge,
                    )
                elif predicate not in self.biolink.predicates:
                    self._add_error("invalid_predicate", f"Invalid predicate '{predicate}'", edge)

            # Knowledge level must be valid per Biolink
            if EdgeModel.knowledge_level.name in edge:
                if edge[EdgeModel.knowledge_level.name] not in self.biolink.knowledge_levels:
                    self._add_error(
                        "invalid_knowledge_level",
                        f"Invalid edge.knowledge_level '{edge[EdgeModel.knowledge_level.name]}'. "
                        f"Valid options are: {self.biolink.knowledge_levels}",
                        edge,
                    )

            # Agent type must be valid per Biolink
            if EdgeModel.agent_type.name in edge:
                if edge[EdgeModel.agent_type.name] not in self.biolink.agent_types:
                    message = (
                        f"Invalid edge.agent_type '{edge[EdgeModel.agent_type.name]}'. "
                        f"Valid options are: {self.biolink.agent_types}"
                    )
                    self._add_error(
                        "invalid_agent_type",
                        message,
                        edge,
                    )

        logging.info(f"Edge validation complete for {edges_path}")

    def _check_file_suffix(self, path: Path) -> None:
        if path.suffix != ".jsonl":
            self._add_error("invalid_file_suffix", f"File must have .jsonl suffix: {path}")

    def _add_error(self, error_type: str, message: str, item: dict | None = None):
        self.summary.add_error(error_type, message, item)

    @staticmethod
    def is_valid_category_format(category: str) -> bool:
        """Check that category is in 'biolink:PascalCase' format."""
        return bool(re.match(r"^biolink:[A-Z][a-zA-Z]*$", category))

    @staticmethod
    def is_valid_predicate_format(predicate: str) -> bool:
        """Check that predicate is in 'biolink:snake_case' format."""
        return bool(re.match(r"^biolink:[a-z][a-z_]*$", predicate))
