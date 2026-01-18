import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import jsonlines

from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.constants import (
    AGENT_TYPE,
    AGGREGATOR_KS,
    ATTRIBUTES,
    CATEGORIES,
    CHEMICAL_FORMULA,
    CONTEXT_QUALIFIER,
    DESCRIPTION,
    EQUIVALENT_IDS,
    EXACT_MASS,
    ID,
    KNOWLEDGE_LEVEL,
    NAME,
    NOT_PROVIDED,
    OBJ_ASPECT_QUALIFIER,
    OBJ_DIRECTION_QUALIFIER,
    OBJECT,
    PREDICATE,
    PRIMARY_KS,
    PROVIDED_BY,
    PUBLICATIONS,
    PUBLICATIONS_INFO,
    QUALIFIED_PREDICATE,
    SUBJECT,
    SUPPORTING_SOURCES,
    SYNONYMS,
    URLS,
)
from kraken.utils.general import clean_text, is_empty, to_list
from kraken.utils.kg_io import (
    stream_edges_from_jsonl,
    stream_edges_from_tsv,
    stream_nodes_from_jsonl,
    stream_nodes_from_tsv,
)


class BaseHarmonizer(ABC):
    """Base class for harmonizing knowledge graph sources into KRAKEN format"""

    @property
    @abstractmethod
    def source_name(self) -> str: ...

    @property
    @abstractmethod
    def source_infores(self) -> str: ...

    list_delimiter: str | None = None

    # Node property name mappings - override when source uses different names
    id_prop: str = ID
    category_prop: str = "category"
    equivalent_ids_prop: str = "xref"
    synonyms_props: set[str] = "synonym"
    url_prop: str = "iri"
    name_prop: str = NAME
    description_prop: str = DESCRIPTION
    chemical_formula_prop: str = CHEMICAL_FORMULA
    exact_mass_prop: str = EXACT_MASS

    # Edge property name mappings - override when source uses different names
    subject_prop: str = SUBJECT
    object_prop: str = OBJECT
    predicate_prop: str = PREDICATE
    primary_ks_prop: str = PRIMARY_KS
    knowledge_level_prop: str = KNOWLEDGE_LEVEL
    agent_type_prop: str = AGENT_TYPE
    qualified_predicate_prop: str = QUALIFIED_PREDICATE
    object_direction_qualifier_prop: str = OBJ_DIRECTION_QUALIFIER
    object_aspect_qualifier_prop: str = OBJ_ASPECT_QUALIFIER
    context_qualifier_prop: str = CONTEXT_QUALIFIER
    supporting_sources_prop: str = SUPPORTING_SOURCES
    publications_prop: str = PUBLICATIONS  # NOTE: this one is also a node prop
    publications_info_prop: str = PUBLICATIONS_INFO

    # Properties to ignore (won't be stored in attributes)
    ignore_node_props: set[str] = set()
    ignore_edge_props: set[str] = set()

    # Rename properties when storing in attributes
    rename_node_attrs: dict[str, str] = {}
    rename_edge_attrs: dict[str, str] = {}

    # Knowledge source overrides
    primary_ks_default_value: str | None = None
    supporting_sources_default_value: str | None = None

    # Properties that should NOT be parsed from delimiter-separated strings (relevant for TSVs only)
    exclude_from_list_parsing: set[str] = set()

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client

        self.core_node_props = {
            self.id_prop,
            self.category_prop,
            self.equivalent_ids_prop,
            self.name_prop,
            self.description_prop,
            self.url_prop,
            self.chemical_formula_prop,
            self.exact_mass_prop,
            self.publications_prop,
        }.union(self.synonyms_props)

        self.core_edge_props = {
            self.subject_prop,
            self.object_prop,
            self.predicate_prop,
            self.primary_ks_prop,
            self.knowledge_level_prop,
            self.agent_type_prop,
            self.qualified_predicate_prop,
            self.object_direction_qualifier_prop,
            self.object_aspect_qualifier_prop,
            self.context_qualifier_prop,
            self.supporting_sources_prop,
            self.publications_prop,
            self.publications_info_prop,
        }

    def harmonize(
        self,
        nodes_output: Path,
        edges_output: Path,
        *,
        input_file: Path | None = None,
        nodes_input: Path | None = None,
        edges_input: Path | None = None,
    ):
        """
        Run full harmonization pipeline.

        Default implementation handles split nodes/edges files.
        Single-file harmonizers should override this method.
        """
        # Validate inputs
        has_separate = nodes_input is not None or edges_input is not None
        has_combined = input_file is not None

        if has_separate and has_combined:
            raise ValueError(f"{self.source_name}: Specify either nodes_input/edges_input OR input_file, not both")
        if has_separate and not (nodes_input and edges_input):
            raise ValueError(f"{self.source_name}: Must specify both nodes_input and edges_input")
        if not has_separate and not has_combined:
            raise ValueError(f"{self.source_name}: Must specify either nodes_input/edges_input or input_file")

        if input_file:
            raise NotImplementedError(
                f"{self.source_name} was called with input_file but does not override harmonize()"
            )

        logging.info(f"Harmonizing {self.source_name}: {nodes_input}, {edges_input} -> {nodes_output}, {edges_output}")
        node_count = self._harmonize_nodes(nodes_input, nodes_output)
        edge_count = self._harmonize_edges(edges_input, edges_output)
        logging.info(f"{self.source_name} harmonization complete: {node_count} nodes, {edge_count} edges")

    def _harmonize_nodes(self, input_path: Path, output_path: Path) -> int:
        count = 0
        with jsonlines.open(output_path, "w") as writer:
            for node in self._stream_nodes(input_path):
                harmonized = self._harmonize_node(node)
                writer.write(harmonized)
                count += 1
        logging.info("Finished harmonizing nodes")
        return count

    def _harmonize_edges(self, input_path: Path, output_path: Path) -> int:
        count = 0
        with jsonlines.open(output_path, "w") as writer:
            for edge in self._stream_edges(input_path):
                harmonized = self._harmonize_edge(edge)
                writer.write(harmonized)
                count += 1
        logging.info("Finished harmonizing edges")
        return count

    def _collect_node_attributes(self, node: dict[str, Any]) -> dict[str, Any]:
        attributes = {}
        for k, v in node.items():
            if k in self.core_node_props or k in self.ignore_node_props or is_empty(v):
                continue
            key = self.rename_node_attrs.get(k, k)
            attributes[key] = v
        return attributes

    def _collect_edge_attributes(self, edge: dict[str, Any]) -> dict[str, Any]:
        attributes = {}
        for k, v in edge.items():
            if k in self.core_edge_props or k in self.ignore_edge_props or is_empty(v):
                continue
            key = self.rename_edge_attrs.get(k, k)
            attributes[key] = v
        return attributes

    def _harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        """Harmonize a single node. Override for source-specific logic."""
        synonyms = set()
        for synonym_prop in self.synonyms_props:
            new_synonyms = to_list(node.get(synonym_prop))
            if new_synonyms:
                synonyms |= set(new_synonyms)

        return self.create_node(
            source_infores=self.source_infores,
            curie=node[self.id_prop],
            categories=self.biolink.filter_to_leaf_categories(node[self.category_prop]),
            provided_by=self.source_infores,
            equivalent_ids=node.get(self.equivalent_ids_prop) if self.equivalent_ids_prop else [node[ID]],
            name=node.get(self.name_prop),
            synonyms=synonyms,
            urls=node.get(self.url_prop),
            description=node.get(self.description_prop),
            chemical_formula=node.get(self.chemical_formula_prop),
            exact_mass=node.get(self.exact_mass_prop),
            publications=node.get(self.publications_prop),
            attributes=self._collect_node_attributes(node),
        )

    def _harmonize_edge(self, edge: dict[str, Any]) -> dict[str, Any]:
        """Harmonize a single edge. Override for source-specific logic."""
        primary_ks = edge[self.primary_ks_prop] if edge.get(self.primary_ks_prop) else self.primary_ks_default_value
        supporting_sources = (
            edge[self.supporting_sources_prop]
            if edge.get(self.supporting_sources_prop)
            else self.supporting_sources_default_value
        )
        return self.create_edge(
            source_infores=self.source_infores,
            subject_id=edge[self.subject_prop],
            object_id=edge[self.object_prop],
            predicate=edge[self.predicate_prop],
            primary_ks=primary_ks,
            knowledge_level=edge.get(self.knowledge_level_prop, NOT_PROVIDED),
            agent_type=edge.get(self.agent_type_prop, NOT_PROVIDED),
            aggregator_ks=self.source_infores,
            qualified_predicate=edge.get(self.qualified_predicate_prop),
            object_direction_qualifier=edge.get(self.object_direction_qualifier_prop),
            object_aspect_qualifier=edge.get(self.object_aspect_qualifier_prop),
            context_qualifier=edge.get(self.context_qualifier_prop),
            supporting_sources=to_list(supporting_sources),
            publications=to_list(edge.get(self.publications_prop, [])),
            publications_info=edge.get(self.publications_info_prop),
            attributes=self._collect_edge_attributes(edge),
        )

    def _stream_nodes(self, input_path: Path | str):
        suffix = Path(input_path).suffix.lower()
        if suffix == ".tsv":
            return stream_nodes_from_tsv(input_path, self.list_delimiter, self.exclude_from_list_parsing)
        elif suffix in (".jsonl", ".jsonlines"):
            return stream_nodes_from_jsonl(input_path)
        else:
            raise ValueError(f"Unknown file format: {suffix}")

    def _stream_edges(self, input_path: Path | str):
        suffix = Path(input_path).suffix.lower()
        if suffix == ".tsv":
            return stream_edges_from_tsv(input_path, self.list_delimiter, self.exclude_from_list_parsing)
        elif suffix in (".jsonl", ".jsonlines"):
            return stream_edges_from_jsonl(input_path)
        else:
            raise ValueError(f"Unknown file format: {suffix}")

    @staticmethod
    def create_node(
        source_infores: str,
        curie: str,
        categories: list[str],
        provided_by: str | list[str],
        equivalent_ids: str | list[str] | None = None,
        name: str | None = None,
        synonyms: list[str] | set[str] | None = None,
        description: str | None = None,
        urls: str | list[str] | None = None,
        chemical_formula: str | list[str] | None = None,
        exact_mass: float | None = None,
        publications: str | list[str] = None,
        attributes: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        if not (curie and categories and provided_by):
            raise ValueError(
                f"Node is missing required field(s): curie={curie}, "
                f"categories={categories}, provided_by={provided_by}"
            )

        # Assemble the node, with properties in a specific order (for convenient review)
        node = {ID: curie}
        if name:
            node[NAME] = clean_text(name)
            # Make sure node's name is in our synonyms list
            if synonyms:
                synonyms = list(set(synonyms) | {name})
            else:
                synonyms = [name]
        node[CATEGORIES] = categories
        if urls:
            node[URLS] = to_list(urls)

        if chemical_formula:
            # Handle case where multiple chemical formulas are given (could be diff nomenclatures, etc.)
            if isinstance(chemical_formula, list):
                if len(chemical_formula) > 1:
                    additional_chemical_formulas = chemical_formula[1:]
                    attributes["additional_chemical_formulas"] = additional_chemical_formulas
                chemical_formula = chemical_formula[0]
            node[CHEMICAL_FORMULA] = chemical_formula

        if exact_mass:
            node[EXACT_MASS] = exact_mass
        if description:
            node[DESCRIPTION] = clean_text(description)

        node[PROVIDED_BY] = to_list(provided_by)  # Convert to list so these will merge

        # Clean up equivalent IDs (remove INCHI, uppercase UNIIs - Molepro has some lowercase)
        cleaned_equiv_ids = set()
        for equiv_id in set(to_list(equivalent_ids) + [curie]):
            equiv_id_upper = equiv_id.upper()
            if not equiv_id_upper.startswith("INCHI:"):
                if equiv_id_upper.startswith("UNII:"):
                    equiv_id = equiv_id_upper
                cleaned_equiv_ids.add(equiv_id)
        node[EQUIVALENT_IDS] = list(cleaned_equiv_ids)

        if synonyms:
            cleaned_synonyms = [clean_text(synonym) for synonym in synonyms if not is_empty(synonym)]
            node[SYNONYMS] = [s for s in cleaned_synonyms if s]

        if publications:
            node[PUBLICATIONS] = to_list(publications)
        if attributes:
            node[ATTRIBUTES] = {source_infores: attributes}

        return node

    @staticmethod
    def create_edge(
        source_infores: str,
        subject_id: str,
        object_id: str,
        predicate: str,
        primary_ks: str | list[str],
        knowledge_level: str,
        agent_type: str | list[str],
        aggregator_ks: str | None = None,
        supporting_sources: list[str] | None = None,
        qualified_predicate: str | None = None,
        object_direction_qualifier: str | None = None,
        object_aspect_qualifier: str | None = None,
        context_qualifier: str | None = None,
        publications: str | list[str] | None = None,
        publications_info: dict[str, Any] | None = None,
        attributes: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        assert subject_id and object_id and predicate and primary_ks and knowledge_level and agent_type

        # Assemble the edge, with properties in a specific order (for convenient review)
        edge = {SUBJECT: subject_id, OBJECT: object_id, PREDICATE: predicate}

        if qualified_predicate:
            edge[QUALIFIED_PREDICATE] = qualified_predicate
        if object_direction_qualifier:
            edge[OBJ_DIRECTION_QUALIFIER] = object_direction_qualifier
        if object_aspect_qualifier:
            edge[OBJ_ASPECT_QUALIFIER] = object_aspect_qualifier
        if context_qualifier:
            edge[CONTEXT_QUALIFIER] = context_qualifier

        # Handle case where multiple primary knowledge sources are given (move others to supporting)
        if isinstance(primary_ks, list):
            if len(primary_ks) > 1:
                supporting_sources = list(set(to_list(supporting_sources) + primary_ks[1:]))
            primary_ks = primary_ks[0]

        # Handle case where multiple agent types are given (just take first, throw others in attributes)
        if isinstance(agent_type, list):
            if len(agent_type) > 1:
                additional_agent_types = agent_type[1:]
                attributes["additional_agent_types"] = additional_agent_types
            agent_type = agent_type[0]

        # Handle improper KLAT
        if knowledge_level == "unspecified":
            knowledge_level = NOT_PROVIDED
        if agent_type == "unspecified":
            agent_type = NOT_PROVIDED

        edge |= {PRIMARY_KS: primary_ks, KNOWLEDGE_LEVEL: knowledge_level, AGENT_TYPE: agent_type}

        if supporting_sources:
            edge[SUPPORTING_SOURCES] = supporting_sources
        if aggregator_ks:
            edge[AGGREGATOR_KS] = to_list(aggregator_ks)  # Convert to list so these will merge
        if publications:
            edge[PUBLICATIONS] = to_list(publications)
        if publications_info:
            edge[PUBLICATIONS_INFO] = publications_info
        if attributes:
            edge[ATTRIBUTES] = {source_infores: attributes}

        return edge

    def validate(self, harmonized_nodes: Path, harmonized_edges: Path):
        logging.info("Validating harmonized nodes/edges...")
        for node in stream_nodes_from_jsonl(harmonized_nodes):
            # TODO: validate nodes!
            pass
        for edge in stream_edges_from_jsonl(harmonized_edges):
            # TODO: validate edges!
            pass
        logging.info("Validation complete! All checks passed.")
