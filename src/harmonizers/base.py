from abc import ABC
from pathlib import Path
from typing import Any
import jsonlines
import logging

from ..utils.constants import *
from ..utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl, stream_nodes_from_tsv, stream_edges_from_tsv
from ..utils.general import to_list, clean_text
from ..utils.biolink_client import BiolinkClient


class BaseHarmonizer(ABC):
    """Base class for harmonizing knowledge graph sources into KRAKEN format"""

    source_name: str
    source_infores: str

    list_delimiter: str = "|"

    # Node property name mappings - override when source uses different names
    id_prop: str = ID
    category_prop: str          # Varies a lot, harmonizer required to set
    equivalent_ids_prop: str    # Varies a lot, harmonizer required to set
    synonyms_props: set[str]    # Varies a lot, harmonizer required to set
    url_prop: str               # Varies a lot, harmonizer required to set
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
    publications_prop: str = PUBLICATIONS
    publications_info_prop: str = PUBLICATIONS_INFO

    # Properties to ignore (won't be stored in attributes)
    ignore_node_props: set[str] = set()
    ignore_edge_props: set[str] = set()

    # Rename properties when storing in attributes
    rename_node_attrs: dict[str, str] = {}
    rename_edge_attrs: dict[str, str] = {}

    empty_values: list[str] = ["", None, []]

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client

        self.core_node_props = {self.id_prop, self.category_prop, self.equivalent_ids_prop,
                                self.name_prop, self.description_prop, self.url_prop,
                                self.chemical_formula_prop, self.exact_mass_prop}.union(self.synonyms_props)
        self.core_edge_props = {self.subject_prop, self.object_prop, self.predicate_prop,
                                self.primary_ks_prop, self.knowledge_level_prop, self.agent_type_prop,
                                self.qualified_predicate_prop, self.object_direction_qualifier_prop,
                                self.object_aspect_qualifier_prop, self.context_qualifier_prop,
                                self.supporting_sources_prop,
                                self.publications_prop, self.publications_info_prop}
        self.could_be_lists = {self.category_prop, self.equivalent_ids_prop,
                               self.supporting_sources_prop, self.publications_prop,
                               self.publications_info_prop}.union(self.synonyms_props)

    def harmonize(
            self,
            nodes_input: Path,
            edges_input: Path,
            nodes_output: Path,
            edges_output: Path,
    ):
        """Run full harmonization pipeline"""
        logging.info(
            f"Harmonizing {self.source_name}: {nodes_input}, {edges_input} -> {nodes_output}, {edges_output}"
        )

        node_count = self._harmonize_nodes(nodes_input, nodes_output)
        edge_count = self._harmonize_edges(edges_input, edges_output)

        logging.info(f"{self.source_name} harmonization complete: {node_count} nodes, {edge_count} edges")

    def _harmonize_nodes(self, input_path: Path, output_path: Path) -> int:
        count = 0
        with jsonlines.open(output_path, 'w') as writer:
            for node in self._stream_nodes(input_path):
                harmonized = self.harmonize_node(node)
                writer.write(harmonized)
                count += 1
        logging.info(f"Finished harmonizing nodes")
        return count

    def _harmonize_edges(self, input_path: Path, output_path: Path) -> int:
        count = 0
        with jsonlines.open(output_path, 'w') as writer:
            for edge in self._stream_edges(input_path):
                harmonized = self.harmonize_edge(edge)
                writer.write(harmonized)
                count += 1
        logging.info(f"Finished harmonizing edges")
        return count

    def collect_node_attributes(self, node: dict[str, Any]) -> dict[str, Any] | None:
        attributes = {}
        for k, v in node.items():
            if k in self.core_node_props or k in self.ignore_node_props or v in self.empty_values:
                continue
            key = self.rename_node_attrs.get(k, k)
            attributes[key] = v
        return {self.source_infores: attributes} if attributes else None

    def collect_edge_attributes(self, edge: dict[str, Any]) -> dict[str, Any] | None:
        attributes = {}
        for k, v in edge.items():
            if k in self.core_edge_props or k in self.ignore_edge_props or v in self.empty_values:
                continue
            key = self.rename_edge_attrs.get(k, k)
            attributes[key] = v
        return {self.source_infores: attributes} if attributes else None

    def harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        """Harmonize a single node. Override for source-specific logic."""
        synonyms = set()
        for synonym_prop in self.synonyms_props:
            new_synonyms = to_list(node.get(synonym_prop))
            if new_synonyms:
                synonyms |= set(new_synonyms)

        return self.create_node(
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
            attributes=self.collect_node_attributes(node)
        )

    def harmonize_edge(self, edge: dict[str, Any]) -> dict[str, Any]:
        """Harmonize a single edge. Override for source-specific logic."""
        return self.create_edge(
            subject_id=edge[self.subject_prop],
            object_id=edge[self.object_prop],
            predicate=edge[self.predicate_prop],
            primary_ks=edge[self.primary_ks_prop],
            knowledge_level=edge.get(self.knowledge_level_prop, NOT_PROVIDED),
            agent_type=edge.get(self.agent_type_prop, NOT_PROVIDED),
            aggregator_ks=self.source_infores,
            qualified_predicate=edge.get(self.qualified_predicate_prop),
            object_direction_qualifier=edge.get(self.object_direction_qualifier_prop),
            object_aspect_qualifier=edge.get(self.object_aspect_qualifier_prop),
            context_qualifier=edge.get(self.context_qualifier_prop),
            supporting_sources=to_list(edge.get(self.supporting_sources_prop)),
            publications=to_list(edge.get(self.publications_prop, [])),
            publications_info=edge.get(self.publications_info_prop),
            attributes=self.collect_edge_attributes(edge)
        )

    def _stream_nodes(self, input_path: Path | str):
        suffix = Path(input_path).suffix.lower()
        if suffix == '.tsv':
            return stream_nodes_from_tsv(input_path, list_delimiter=self.list_delimiter, could_be_list=self.could_be_lists)
        elif suffix in ('.jsonl', '.jsonlines'):
            return stream_nodes_from_jsonl(input_path)
        else:
            raise ValueError(f"Unknown file format: {suffix}")

    def _stream_edges(self, input_path: Path | str):
        suffix = Path(input_path).suffix.lower()
        if suffix == '.tsv':
            return stream_edges_from_tsv(input_path, list_delimiter=self.list_delimiter, could_be_list=self.could_be_lists)
        elif suffix in ('.jsonl', '.jsonlines'):
            return stream_edges_from_jsonl(input_path)
        else:
            raise ValueError(f"Unknown file format: {suffix}")


    @staticmethod
    def create_node(curie: str,
                    categories: list[str],
                    provided_by: str | list[str],
                    equivalent_ids: str | list[str] | None = None,
                    name: str | None = None,
                    synonyms: list[str] | set[str] | None = None,
                    description: str | None = None,
                    urls: str | list[str] | None = None,
                    chemical_formula: str | None = None,
                    exact_mass: float | None = None,
                    attributes: dict[str, Any] | None = None) -> dict[str, Any]:
        if not (curie and categories and provided_by):
            raise ValueError(f"Node is missing required field(s): curie={curie}, "
                             f"categories={categories}, provided_by={provided_by}")

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
            node[CHEMICAL_FORMULA] = chemical_formula
        if exact_mass:
            node[EXACT_MASS] = exact_mass
        if description:
            node[DESCRIPTION] = clean_text(description)

        node[PROVIDED_BY] = to_list(provided_by)  # Convert to list so these will merge
        node[EQUIVALENT_IDS] = list(set(to_list(equivalent_ids) + [curie]))

        if synonyms:
            cleaned_synonyms = [clean_text(synonym) for synonym in synonyms if synonym]
            node[SYNONYMS] = [s for s in cleaned_synonyms if s]

        if attributes:
            node[ATTRIBUTES] = attributes

        return node

    @staticmethod
    def create_edge(subject_id: str,
                    object_id: str,
                    predicate: str,
                    primary_ks: str,
                    knowledge_level: str,
                    agent_type: str,
                    aggregator_ks: str | None = None,
                    supporting_sources: list[str] | None = None,
                    qualified_predicate: str | None = None,
                    object_direction_qualifier: str | None = None,
                    object_aspect_qualifier: str | None = None,
                    context_qualifier: str | None = None,
                    publications: list[str] | None = None,
                    publications_info: dict[str, Any] | None = None,
                    attributes: dict[str, Any] | None = None) -> dict[str, Any]:
        assert subject_id and object_id and predicate and primary_ks and knowledge_level and agent_type

        # Assemble the edge, with properties in a specific order (for convenient review)
        edge = {SUBJECT: subject_id,
                OBJECT: object_id,
                PREDICATE: predicate}

        if qualified_predicate:
            edge[QUALIFIED_PREDICATE] = qualified_predicate
        if object_direction_qualifier:
            edge[OBJ_DIRECTION_QUALIFIER] = object_direction_qualifier
        if object_aspect_qualifier:
            edge[OBJ_ASPECT_QUALIFIER] = object_aspect_qualifier
        if context_qualifier:
            edge[CONTEXT_QUALIFIER] = context_qualifier

        edge |= {PRIMARY_KS: primary_ks,
                 KNOWLEDGE_LEVEL: knowledge_level,
                 AGENT_TYPE: agent_type}

        if supporting_sources:
            edge[SUPPORTING_SOURCES] = supporting_sources
        if aggregator_ks:
            edge[AGGREGATOR_KS] = to_list(aggregator_ks)  # Convert to list so these will merge
        if publications:
            edge[PUBLICATIONS] = publications
        if publications_info:
            edge[PUBLICATIONS_INFO] = publications_info
        if attributes:
            edge[ATTRIBUTES] = attributes

        return edge