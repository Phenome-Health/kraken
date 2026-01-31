"""Validator for KRAKEN-ified JSONL files."""

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from kraken.biolink_client import BiolinkClient
from kraken.schema import EdgeModel, NodeModel, PropertyDef
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl


@dataclass
class ValidationError:
    """A single validation error."""

    message: str
    item: dict | None = None


@dataclass
class ValidationSummary:
    """Collects and summarizes validation errors."""

    errors_by_subtype: dict[str, dict[str | None, list[ValidationError]]] = field(
        default_factory=lambda: defaultdict(lambda: defaultdict(list))
    )
    max_examples_per_type: int = 3
    max_examples_per_subtype: int = 1

    def add_error(self, error_type: str, message: str, item: dict | None = None, subtype: str | None = None):
        self.errors_by_subtype[error_type][subtype].append(ValidationError(message, item))

    @property
    def total_errors(self) -> int:
        return sum(len(errors) for subtypes in self.errors_by_subtype.values() for errors in subtypes.values())

    @property
    def has_errors(self) -> bool:
        return self.total_errors > 0

    def get_summary(self) -> str:
        if not self.has_errors:
            return "No validation errors found."

        lines = [
            f"Validation failed with {self.total_errors} total error(s) across {len(self.errors_by_subtype)} type(s):\n"
        ]

        for error_type, subtypes in sorted(self.errors_by_subtype.items()):
            type_total = sum(len(errors) for errors in subtypes.values())
            lines.append(f"[{error_type}] - {type_total} occurrence(s)")

            # Check if we have real subtypes (not just None)
            has_subtypes = any(subtype is not None for subtype in subtypes.keys())

            if has_subtypes:
                for subtype, subtype_errors in sorted(
                    ((k, v) for k, v in subtypes.items() if k is not None), key=lambda x: -len(x[1])
                ):
                    lines.append(f"\n    • {subtype}: {len(subtype_errors)}")
                    for error in subtype_errors[: self.max_examples_per_subtype]:
                        if error.item:
                            lines.append(f"      Example: {error.item}")
            else:
                # No subtypes, show examples directly
                all_errors = subtypes[None]
                for error in all_errors[: self.max_examples_per_type]:
                    lines.append(f"    • {error.message}")
                    if error.item:
                        lines.append(f"      Item: {error.item}")
                if len(all_errors) > self.max_examples_per_type:
                    lines.append(f"    ... and {len(all_errors) - self.max_examples_per_type} more")

        return "\n".join(lines)


class KrakenValidator:
    """Validates KRAKEN-harmonized JSONL node and edge files."""

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client
        self.summary = ValidationSummary()
        self.node_ids = set()

    def validate(
        self, nodes_path: Path, edges_path: Path, source_infores: str | None = None, integrated: bool = False
    ) -> None:
        """
        Validate KRAKEN-harmonized node and edge files.

        Raises:
            ValueError: If any validation check fails.
        """
        logging.info(f"Validating {nodes_path} and {edges_path}...")

        # Reset summary for fresh validation run
        self.summary = ValidationSummary()

        self.validate_nodes(nodes_path, source_infores, integrated)
        self.validate_edges(edges_path, integrated)

        # TODO: attributes (nodes and edges), types within lists, non-core edge props, aggregator ks.. qualifiers?
        # TODO: name, synonyms, none strings?

        if self.summary.has_errors:
            raise ValueError(self.summary.get_summary())

        logging.info("Validation complete! All checks passed. ✅")

    def validate_nodes(self, nodes_path: Path, source_infores: str | None = None, integrated: bool = False) -> None:
        logging.info(f"Starting to validate {nodes_path}")

        self._check_file_suffix(nodes_path)

        required_props = [p for p in NodeModel.all_properties().values() if p.required]
        all_props = {p.name for p in NodeModel.all_properties().values()}
        merged_node_count = 0

        for node in stream_nodes_from_jsonl(nodes_path):
            # Record this node ID for use during edge validation
            if node.get(NodeModel.id.name):
                self.node_ids.add(node[NodeModel.id.name])

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
                    if not isinstance(value, prop.dtype):
                        self._add_error(
                            "invalid_type",
                            f"Field node.{prop.name} has type {type(value).__name__}, expected {prop.dtype.__name__}",
                            node,
                        )
                    else:
                        self._validate_inner_types(prop, value, node, "node")

            # Check for unexpected properties
            for prop_name in node.keys():
                if prop_name not in all_props:
                    self._add_error(
                        "unexpected_property",
                        f"Unexpected property 'node.{prop_name}' not defined in NodeModel",
                        node,
                        subtype=prop_name,
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
                            subtype=category,
                        )
                    elif not (category in self.biolink.categories or self.biolink.toolkit.is_mixin(category)):
                        # NOTE: biolink has bug where some category mixins are not descendants of NamedThing
                        self._add_error(
                            "invalid_category",
                            f"Category '{category}' does not exist in Biolink",
                            node,
                            subtype=category,
                        )

            if integrated:
                if len(node[NodeModel.provided_by.name]) > 1:
                    merged_node_count += 1
                    # Print out the first few merged nodes
                    if merged_node_count < 3:
                        logging.info(f"Merged node example: {node}")

        # Ensure there are some merged nodes if these are integrated kraken files
        if integrated:
            if not merged_node_count:
                message = "No nodes merged from multiple kraken sources detected in integrated files"
                self._add_error("no_merged_nodes", message)
            else:
                logging.info(f"Detected {merged_node_count} nodes merged from multiple kraken sources")

    def validate_edges(self, edges_path: Path, integrated: bool = False) -> None:
        logging.info(f"Starting to validate {edges_path}")

        self._check_file_suffix(edges_path)

        required_props = [p for p in EdgeModel.all_properties().values() if p.required]
        all_props = {p.name for p in EdgeModel.all_properties().values()}
        merged_edge_count = True

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
                    if not isinstance(value, prop.dtype):
                        self._add_error(
                            "invalid_type",
                            f"Field edge.{prop.name} has type {type(value).__name__}, expected {prop.dtype.__name__}",
                            edge,
                        )
                    else:
                        self._validate_inner_types(prop, value, edge, "edge")

            # Orphan edge check
            subject_id = edge.get(EdgeModel.subject.name)
            object_id = edge.get(EdgeModel.object.name)
            if subject_id and subject_id not in self.node_ids:
                self._add_error(
                    "orphan_edge",
                    f"Edge subject '{subject_id}' does not exist in nodes",
                    edge,
                    subtype=EdgeModel.subject.name,
                )
            if object_id and object_id not in self.node_ids:
                self._add_error(
                    "orphan_edge",
                    f"Edge object '{object_id}' does not exist in nodes",
                    edge,
                    subtype=EdgeModel.object.name,
                )

            # Check for unexpected properties
            for prop_name in edge.keys():
                if prop_name not in all_props:
                    self._add_error(
                        "unexpected_property",
                        f"Unexpected property 'edge.{prop_name}' not defined in EdgeModel",
                        edge,
                        subtype=prop_name,
                    )

            # Predicate must be valid Biolink predicate
            if EdgeModel.predicate.name in edge:
                predicate = edge[EdgeModel.predicate.name]
                if not self.is_valid_predicate_format(predicate):
                    self._add_error(
                        "invalid_predicate_format",
                        f"Predicate '{predicate}' is not in proper biolink:snake_case format",
                        edge,
                        subtype=predicate,
                    )
                elif predicate not in self.biolink.predicates:
                    self._add_error(
                        "invalid_predicate",
                        f"Invalid predicate '{predicate}'",
                        edge,
                        subtype=predicate,
                    )

            # Knowledge level must be valid per Biolink
            if EdgeModel.knowledge_level.name in edge:
                knowledge_level = edge[EdgeModel.knowledge_level.name]
                if knowledge_level not in self.biolink.knowledge_levels:
                    self._add_error(
                        "invalid_knowledge_level",
                        f"Invalid edge.knowledge_level '{knowledge_level}'. "
                        f"Valid options are: {self.biolink.knowledge_levels}",
                        edge,
                        subtype=knowledge_level,
                    )

            # Agent type must be valid per Biolink
            if EdgeModel.agent_type.name in edge:
                agent_type = edge[EdgeModel.agent_type.name]
                if agent_type not in self.biolink.agent_types:
                    primary_ks = edge.get(EdgeModel.primary_ks.name, "unknown")
                    message = (
                        f"Invalid edge.agent_type '{agent_type}'. " f"Valid options are: {self.biolink.agent_types}"
                    )
                    self._add_error(
                        "invalid_agent_type",
                        message,
                        edge,
                        subtype=f"{agent_type} (source: {primary_ks})",
                    )

            # Record whether this is a merged edge from multiple aggregators
            if EdgeModel.aggregator_ks.name in edge:
                if len(edge[EdgeModel.aggregator_ks.name]) > 1:
                    merged_edge_count += 1
                    # Print out the first few merged edges
                    if merged_edge_count < 3:
                        logging.info(f"Merged edge example: {edge}")

        # Ensure there are some merged edges if these are integrated kraken files
        if integrated:
            if not merged_edge_count:
                message = "No edges merged from multiple aggregators detected in integrated files"
                self._add_error("no_merged_edges", message)
            else:
                logging.info(f"Detected {merged_edge_count} edges merged from multiple aggregators")

        logging.info(f"Edge validation complete for {edges_path}")

    def _check_file_suffix(self, path: Path) -> None:
        if path.suffix != ".jsonl":
            self._add_error("invalid_file_suffix", f"File must have .jsonl suffix: {path}")

    def _add_error(self, error_type: str, message: str, item: dict | None = None, subtype: str | None = None):
        self.summary.add_error(error_type, message, item, subtype)

    @staticmethod
    def is_valid_category_format(category: str) -> bool:
        """Check that category is in 'biolink:PascalCase' format."""
        return bool(re.match(r"^biolink:[A-Z][a-zA-Z]*$", category))

    @staticmethod
    def is_valid_predicate_format(predicate: str) -> bool:
        """Check that predicate is in 'biolink:snake_case' format."""
        return bool(re.match(r"^biolink:[a-z][a-z_]*$", predicate))

    def _validate_inner_types(self, prop: PropertyDef, value: any, item: dict, item_label: str) -> None:
        """Validate types of items within a list."""
        if prop.inner_type is None or not isinstance(value, list):
            return

        for i, inner_value in enumerate(value):
            if not isinstance(inner_value, prop.inner_type):
                self._add_error(
                    "invalid_inner_type",
                    f"Field {item_label}.{prop.name}[{i}] has type {type(inner_value).__name__}, "
                    f"expected {prop.inner_type.__name__}",
                    item,
                )
