import logging
from abc import ABC, abstractmethod
from collections import defaultdict
from pathlib import Path
from typing import Any

import jsonlines
from biomapper2.core.normalizer import Normalizer

from kraken.biolink_client import BiolinkClient
from kraken.utils.constants import (
    EDGE_AGENT_TYPE,
    EDGE_AGGREGATOR_KS,
    EDGE_ATTRIBUTES,
    EDGE_KNOWLEDGE_LEVEL,
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_PRIMARY_KS,
    EDGE_PUBLICATIONS,
    EDGE_PUBLICATIONS_INFO,
    EDGE_QUALIFIERS,
    EDGE_SUBJECT,
    EDGE_SUPPORTING_SOURCES,
    INFORES_PREFIX,
    NODE_ATTRIBUTES,
    NODE_CATEGORIES,
    NODE_CHEMICAL_FORMULA,
    NODE_DESCRIPTION,
    NODE_EQUIVALENT_IDS,
    NODE_EXACT_MASS,
    NODE_ID,
    NODE_NAME,
    NODE_PROVIDED_BY,
    NODE_PUBLICATIONS,
    NODE_SYNONYMS,
    NODE_TAXA,
    NODE_URLS,
    NOT_PROVIDED,
    UNRELIABLE_PUBLICATION_PRIMARY_KS,
)
from kraken.utils.general import clean_text, is_empty, to_list
from kraken.utils.kg_io import (
    fix_repeated_prefix,
    split_curie,
    stream_edges_from_jsonl,
    stream_edges_from_tsv,
    stream_nodes_from_jsonl,
    stream_nodes_from_tsv,
)

# TRAPI-style edge provenance: an edge carrying this field lists its knowledge sources as objects of
# {resource_id, resource_role}, which we parse instead of the flat primary_ks/supporting_sources props.
TRAPI_SOURCES_FIELD = "sources"
TRAPI_SOURCE_ROLES = {"primary_knowledge_source", "aggregator_knowledge_source", "supporting_data_source"}


class BaseHarmonizer(ABC):
    """Base class for harmonizing knowledge graph sources into KRAKEN format"""

    @property
    @abstractmethod
    def source_infores(self) -> str: ...

    @property
    def source_name(self) -> str:
        """Human-readable label used in log/error messages, derived from source_infores by stripping any
        'infores:' prefix (e.g. 'infores:rtx-kg2' -> 'rtx-kg2'; a bare id like 'nih-cde' is used as-is)."""
        return self.source_infores.removeprefix(f"{INFORES_PREFIX}:")

    list_delimiter: str | None = None  # Relevant only for TSVs/CSVs
    is_aggregator: bool = False

    # Node property name mappings - override when source uses different names
    id_prop: str = NODE_ID
    category_prop: str = "category"
    equivalent_ids_prop: str = "xref"
    synonyms_props: set[str] = {"synonym"}
    taxon_props: set[str] = set()  # source field(s) to pull taxon CURIE(s) from (unioned), e.g. {"taxon", "in_taxon"}
    url_prop: str = "iri"
    name_prop: str = NODE_NAME
    description_prop: str = NODE_DESCRIPTION
    chemical_formula_prop: str = NODE_CHEMICAL_FORMULA
    exact_mass_prop: str = NODE_EXACT_MASS

    # Edge property name mappings - override when source uses different names
    subject_prop: str = EDGE_SUBJECT
    object_prop: str = EDGE_OBJECT
    predicate_prop: str = EDGE_PREDICATE
    primary_ks_prop: str = EDGE_PRIMARY_KS
    knowledge_level_prop: str = EDGE_KNOWLEDGE_LEVEL
    agent_type_prop: str = EDGE_AGENT_TYPE
    supporting_sources_prop: str = EDGE_SUPPORTING_SOURCES
    publications_prop: str = EDGE_PUBLICATIONS  # NOTE: this one is also a node prop
    publications_info_prop: str = EDGE_PUBLICATIONS_INFO

    # Properties to ignore (won't be stored in attributes)
    # Source field names to drop entirely from harmonized output -- neither mapped to a top-level property nor
    # retained in attributes (applied before any prop mapping). Don't list required fields (id/category for nodes;
    # subject/object/predicate/primary_ks for edges) -- harmonization errors out if a required field is dropped.
    ignore_node_props: set[str] = set()
    ignore_edge_props: set[str] = set()

    # Rename properties when storing in attributes/qualifiers
    rename_node_attrs: dict[str, str] = {}
    rename_edge_attrs_or_quals: dict[str, str] = {}

    # Knowledge source defaults
    primary_ks_default_value: str | None = None
    supporting_sources_default_value: str | None = None

    # Overrides for specific categories, predicates, or agent types
    predicate_overrides: dict[str, str] = dict()
    category_overrides: dict[str, str] = dict()
    agent_type_overrides: dict[str, dict[str, str]] = dict()  # Organized by primary KS

    # Primary knowledge sources to skip (these edges will NOT be included)
    primary_ks_exclusions: set = set()

    # Drop edges asserting a negation (negated=true). KRAKEN has no way to represent negation, so ingesting
    # such an edge would wrongly read as a positive assertion. Applies to all sources; set False to keep them.
    drop_negated_edges: bool = True
    negated_prop: str = "negated"

    # Properties that should NOT be parsed from delimiter-separated strings (relevant for TSVs only)
    exclude_from_list_parsing: set[str] = set()

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client
        # Set up biomapper2's normalizer, so we can normalize curies as needed
        self.normalizer = Normalizer(biolink_version=self.biolink.version)
        self.unrecognized_vocabs = set()
        self.prefixes_with_invalid_ids = defaultdict(int)
        self.invalid_curies = set()
        self.normalized_id_map = dict()
        self.unrecognized_source_roles = set()
        self.stripped_publications_count = 0

        self.core_node_props = (
            {
                self.id_prop,
                self.category_prop,
                self.equivalent_ids_prop,
                self.name_prop,
                self.description_prop,
                self.url_prop,
                self.chemical_formula_prop,
                self.exact_mass_prop,
                self.publications_prop,
            }
            .union(self.synonyms_props)
            .union(self.taxon_props)
        )

        self.core_edge_props = {
            self.subject_prop,
            self.object_prop,
            self.predicate_prop,
            self.primary_ks_prop,
            self.knowledge_level_prop,
            self.agent_type_prop,
            self.supporting_sources_prop,
            self.publications_prop,
            self.publications_info_prop,
        }
        # NOTE: the TRAPI `sources` field is deliberately NOT in core_edge_props, so the full raw list is
        # retained in the edge's attributes (in addition to being flattened into primary/aggregator/supporting).

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

        if self.invalid_curies:
            logging.warning(
                f"A total of {len(self.invalid_curies)} nodes had IDs that are not curies (left them as they are). "
                f"First 50 are: {list(self.invalid_curies)[:10]}"
            )
        if self.unrecognized_vocabs:
            logging.warning(
                f"biomapper2 failed to recognize {len(self.unrecognized_vocabs)} vocabs: {self.unrecognized_vocabs}"
            )
        if self.prefixes_with_invalid_ids:
            logging.warning(
                f"some IDs failed validation in biomapper2 (left them as they are) - counts by prefix are: "
                f"{dict(sorted(self.prefixes_with_invalid_ids.items(), key=lambda x: x[1], reverse=True))}"
            )

        return count

    def _harmonize_edges(self, input_path: Path, output_path: Path) -> int:
        count = 0
        excluded_count = 0
        negated_count = 0
        self.stripped_publications_count = 0
        with jsonlines.open(output_path, "w") as writer:
            for edge in self._stream_edges(input_path):
                # Skip negated edges (KRAKEN can't represent negation, so they'd read as positive assertions)
                if self.drop_negated_edges and edge.get(self.negated_prop):
                    negated_count += 1
                    continue
                # Skip edges from primary knowledge sources marked for exclusion
                primary_kses = set(to_list(edge.get(self.primary_ks_prop)))
                if self.primary_ks_exclusions and primary_kses.issubset(self.primary_ks_exclusions):
                    excluded_count += 1
                else:
                    # Add this edge
                    harmonized = self._harmonize_edge(edge)
                    writer.write(harmonized)
                    count += 1
        logging.info("Finished harmonizing edges.")
        if negated_count:
            logging.info(f"Dropped {negated_count} negated edges (negated=true).")
        if self.stripped_publications_count:
            logging.info(
                f"Dropped publications from {self.stripped_publications_count} edges whose primary knowledge "
                f"source has unreliable publications ({sorted(UNRELIABLE_PUBLICATION_PRIMARY_KS)})."
            )
        if excluded_count:
            logging.info(
                f"Excluded {excluded_count} edges that came from these primary "
                f"knowledge sources: {self.primary_ks_exclusions}."
            )
        return count

    def _collect_node_attributes(self, node: dict[str, Any]) -> dict[str, Any]:
        attributes = {}
        for k, v in node.items():
            if k in self.core_node_props or k in self.ignore_node_props or is_empty(v):
                continue
            key = self.rename_node_attrs.get(k, k)
            attributes[key] = v
        return attributes

    def _collect_edge_attributes_and_qualifiers(self, edge: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
        attributes, qualifiers = {}, {}
        for k, v in edge.items():
            if k in self.core_edge_props or k in self.ignore_edge_props or is_empty(v):
                continue
            key = self.rename_edge_attrs_or_quals.get(k, k)
            if "qualifie" in key.lower():
                qualifiers[key.lower()] = v
            else:
                attributes[key] = v
        return attributes, qualifiers

    def _harmonize_node(self, node: dict[str, Any]) -> dict[str, Any]:
        """Harmonize a single node. Override for source-specific logic."""
        # Drop ignored source fields up front so they reach neither a top-level property nor attributes
        node = {k: v for k, v in node.items() if k not in self.ignore_node_props}

        # Grab all synonyms as applicable
        synonyms = set()
        for synonym_prop in self.synonyms_props:
            new_synonyms = to_list(node.get(synonym_prop))
            if new_synonyms:
                synonyms |= set(new_synonyms)

        # Collect taxon CURIE(s) from any configured taxon field(s), unioned
        taxa = set()
        for taxon_prop in self.taxon_props:
            taxa |= set(to_list(node.get(taxon_prop)))

        # TODO: is it worth having a distinction between these and create_node/edge functions?

        return self.create_node(
            curie=node[self.id_prop],
            categories=node[self.category_prop],
            provided_by=self.source_infores,
            equivalent_ids=node.get(self.equivalent_ids_prop) if self.equivalent_ids_prop else [node[NODE_ID]],
            name=node.get(self.name_prop),
            synonyms=synonyms,
            urls=node.get(self.url_prop),
            description=node.get(self.description_prop),
            chemical_formula=node.get(self.chemical_formula_prop),
            exact_mass=node.get(self.exact_mass_prop),
            publications=node.get(self.publications_prop),
            taxa=taxa,
            attributes=self._collect_node_attributes(node),
        )

    def _parse_trapi_sources(self, edge: dict[str, Any]) -> tuple[str, list[str], list[str]]:
        """Extract (primary_ks, aggregator_kses, supporting_sources) from an edge's TRAPI-style `sources`
        list (each entry a {resource_id, resource_role} object). Raises on anything that isn't valid TRAPI,
        so malformed provenance halts the build rather than being silently mis-recorded. upstream_resource_ids
        (the provenance chain) aren't representable in KRAKEN's flat model and are dropped here -- but the full
        raw `sources` list is retained in the edge's attributes."""
        sources = edge[TRAPI_SOURCES_FIELD]
        edge_id = edge.get("id", "?")
        if not isinstance(sources, list):
            raise ValueError(
                f"{self.source_name}: edge {edge_id!r} has a non-list '{TRAPI_SOURCES_FIELD}': {sources!r}"
            )
        primary_ks, aggregator_kses, supporting_sources = None, [], []
        for source in sources:
            if not isinstance(source, dict) or not source.get("resource_id") or not source.get("resource_role"):
                raise ValueError(
                    f"{self.source_name}: edge {edge_id!r} has a malformed '{TRAPI_SOURCES_FIELD}' entry "
                    f"(each needs resource_id + resource_role): {source!r}"
                )
            resource_id, role = source["resource_id"], source["resource_role"]
            if role not in TRAPI_SOURCE_ROLES:
                # Not one of our three captured roles: skip flattening it (the full source entry is still
                # retained in the edge's attributes). Warn once per unrecognized role to avoid log spam.
                if role not in self.unrecognized_source_roles:
                    self.unrecognized_source_roles.add(role)
                    logging.warning(
                        f"{self.source_name}: unrecognized resource_role {role!r} in '{TRAPI_SOURCES_FIELD}' "
                        f"(not one of {sorted(TRAPI_SOURCE_ROLES)}); it won't be flattened into "
                        f"primary/aggregator/supporting, but the raw source is retained in attributes."
                    )
                continue
            if role == "primary_knowledge_source":
                if primary_ks is not None:
                    raise ValueError(f"{self.source_name}: edge {edge_id!r} has multiple primary_knowledge_sources")
                primary_ks = resource_id
            elif role == "aggregator_knowledge_source":
                aggregator_kses.append(resource_id)
            else:  # supporting_data_source
                supporting_sources.append(resource_id)
        if primary_ks is None:
            raise ValueError(
                f"{self.source_name}: edge {edge_id!r} has no primary_knowledge_source in '{TRAPI_SOURCES_FIELD}'"
            )
        return primary_ks, aggregator_kses, supporting_sources

    def _harmonize_edge(self, edge: dict[str, Any]) -> dict[str, Any]:
        """Harmonize a single edge. Override for source-specific logic."""
        # Drop ignored source fields up front so they reach neither a top-level property nor attributes
        edge = {k: v for k, v in edge.items() if k not in self.ignore_edge_props}

        # Determine knowledge sources from a TRAPI-style `sources` list when the edge has one, else flat props.
        # (The raw `sources` list is also retained in the edge's attributes -- see _collect_edge_attributes_*.)
        if edge.get(TRAPI_SOURCES_FIELD):
            primary_ks, aggregator_ks, supporting_sources = self._parse_trapi_sources(edge)
            if self.is_aggregator:  # also record this KG itself as an aggregator of the edge
                aggregator_ks = list(dict.fromkeys(aggregator_ks + [self.source_infores]))
        else:
            primary_ks = edge[self.primary_ks_prop] if edge.get(self.primary_ks_prop) else self.primary_ks_default_value
            supporting_sources = (
                edge[self.supporting_sources_prop]
                if edge.get(self.supporting_sources_prop)
                else self.supporting_sources_default_value
            )
            aggregator_ks = self.source_infores if self.is_aggregator else None
        attributes, qualifiers = self._collect_edge_attributes_and_qualifiers(edge)

        # Drop publications from sources whose publication lists are unreliable (see the constant's docstring)
        primary_ks_id = primary_ks[0] if isinstance(primary_ks, list) else primary_ks
        if primary_ks_id in UNRELIABLE_PUBLICATION_PRIMARY_KS:
            publications, publications_info = [], None
            self.stripped_publications_count += 1
        else:
            publications = to_list(edge.get(self.publications_prop, []))
            publications_info = edge.get(self.publications_info_prop)

        return self.create_edge(
            subject_id=edge[self.subject_prop],
            object_id=edge[self.object_prop],
            predicate=edge[self.predicate_prop],
            primary_ks=primary_ks,
            knowledge_level=edge.get(self.knowledge_level_prop, NOT_PROVIDED),
            agent_type=edge.get(self.agent_type_prop, NOT_PROVIDED),
            aggregator_ks=aggregator_ks,
            supporting_sources=to_list(supporting_sources),
            publications=publications,
            publications_info=publications_info,
            qualifiers=qualifiers,
            attributes=attributes,
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

    def create_node(
        self,
        curie: str,
        categories: str | list[str],
        provided_by: str | list[str],
        equivalent_ids: str | list[str] | None = None,
        name: str | None = None,
        synonyms: list[str] | set[str] | None = None,
        description: str | None = None,
        urls: str | list[str] | None = None,
        chemical_formula: str | list[str] | None = None,
        exact_mass: float | None = None,
        publications: str | list[str] = None,
        taxa: str | list[str] | set[str] | None = None,
        attributes: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        if not (curie and categories and provided_by):
            raise ValueError(
                f"Node is missing required field(s): curie={curie}, "
                f"categories={categories}, provided_by={provided_by}"
            )
        if attributes is None:
            attributes = {}

        # Assemble the node, with properties in a specific order (for convenient review)
        curie_normalized = self.normalize_curie(curie)
        node = {NODE_ID: curie_normalized}
        if name:
            node[NODE_NAME] = clean_text(name)
            # Make sure node's name is in our synonyms list
            if synonyms:
                synonyms = list(set(synonyms) | {name})
            else:
                synonyms = [name]

        # Clean categories and override as applicable
        categories_cleaned = [fix_repeated_prefix(category) for category in to_list(categories)]
        categories = [self.category_overrides.get(category, category) for category in categories_cleaned]
        node[NODE_CATEGORIES] = self.biolink.filter_to_leaf_categories(categories)

        if urls:
            node[NODE_URLS] = to_list(urls)

        if chemical_formula:
            # Handle case where multiple chemical formulas are given (could be diff nomenclatures, etc.)
            if isinstance(chemical_formula, list):
                if len(chemical_formula) > 1:
                    additional_chemical_formulas = chemical_formula[1:]
                    attributes["additional_chemical_formulas"] = additional_chemical_formulas
                chemical_formula = chemical_formula[0]
            node[NODE_CHEMICAL_FORMULA] = chemical_formula

        if exact_mass:
            node[NODE_EXACT_MASS] = float(exact_mass)
        if description:
            node[NODE_DESCRIPTION] = clean_text(description)

        node[NODE_PROVIDED_BY] = to_list(provided_by)  # Convert to list so these will merge

        # Clean up equivalent IDs (normalize and remove INCHI IDs - big and not very helpful)
        cleaned_equiv_ids = set()
        for equiv_id in to_list(equivalent_ids):
            if not equiv_id.startswith("INCHI:"):
                normalized_equiv_id = self.normalize_curie(equiv_id)
                cleaned_equiv_ids.add(normalized_equiv_id)
        equivalent_ids_final = list(cleaned_equiv_ids | {curie_normalized})
        node[NODE_EQUIVALENT_IDS] = equivalent_ids_final

        if taxa:
            node[NODE_TAXA] = list(dict.fromkeys(to_list(taxa)))  # list so taxa union-merge across sources

        if synonyms:
            cleaned_synonyms = [clean_text(synonym) for synonym in synonyms if not is_empty(synonym)]
            node[NODE_SYNONYMS] = [s for s in cleaned_synonyms if s]

        if publications:
            node[NODE_PUBLICATIONS] = to_list(publications)
        if attributes:
            node[NODE_ATTRIBUTES] = {self.source_infores: attributes}

        return node

    def create_edge(
        self,
        subject_id: str,
        object_id: str,
        predicate: str,
        primary_ks: str | list[str],
        knowledge_level: str,
        agent_type: str | list[str],
        aggregator_ks: str | list[str] | None = None,
        supporting_sources: list[str] | None = None,
        publications: str | list[str] | None = None,
        publications_info: dict[str, Any] | None = None,
        qualifiers: dict[str, Any] | None = None,
        attributes: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        assert subject_id and object_id and predicate and primary_ks and knowledge_level and agent_type
        if attributes is None:
            attributes = {}

        # Clean predicate and override as applicable
        predicate_cleaned = fix_repeated_prefix(predicate)
        predicate = self.predicate_overrides.get(predicate_cleaned, predicate_cleaned)

        # Assemble the edge, with properties in a specific order (for convenient review)
        edge = {
            EDGE_SUBJECT: self.normalize_curie(subject_id),
            EDGE_OBJECT: self.normalize_curie(object_id),
            EDGE_PREDICATE: predicate,
        }

        if qualifiers:
            edge[EDGE_QUALIFIERS] = qualifiers

        # Handle case where multiple primary knowledge sources are given (move others to supporting)
        if isinstance(primary_ks, list):
            if len(primary_ks) > 1:
                supporting_sources = list(set(to_list(supporting_sources) + primary_ks[1:]))
            primary_ks = primary_ks[0]

        # Override agent type(s) as applicable
        if self.agent_type_overrides:
            raw_agent_types = to_list(agent_type)
            agent_types = [
                self.agent_type_overrides.get(raw_agent_type, dict()).get(primary_ks, raw_agent_type)
                for raw_agent_type in raw_agent_types
            ]
            agent_type = list(dict.fromkeys(agent_types))  # Deduplicate while preserving order

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

        edge |= {EDGE_PRIMARY_KS: primary_ks, EDGE_KNOWLEDGE_LEVEL: knowledge_level, EDGE_AGENT_TYPE: agent_type}

        if supporting_sources:
            edge[EDGE_SUPPORTING_SOURCES] = supporting_sources
        if aggregator_ks:
            edge[EDGE_AGGREGATOR_KS] = to_list(aggregator_ks)  # Convert to list so these will merge
        if publications:
            edge[EDGE_PUBLICATIONS] = to_list(publications)
        if publications_info:
            edge[EDGE_PUBLICATIONS_INFO] = publications_info

        if attributes:
            edge[EDGE_ATTRIBUTES] = {self.source_infores: attributes}

        return edge

    def normalize_curie(self, curie: str) -> str:
        # TODO: Eventually run all curies through biomapper, but some bug fixes are needed first
        # TODO: Temporarily we'll just run it on molepro's known problem curies
        if ":" in curie and curie.split(":")[0].upper() in {"CHEMBL.COMPOUND", "CHEMBL.TARGET", "UNII", "KEGG"}:
            # Returned the cached mapping if we've seen this curie before
            if curie in self.normalized_id_map:
                return self.normalized_id_map[curie]

            try:
                prefix, local_id = split_curie(curie)
            except Exception:
                self.invalid_curies.add(curie)
                self.normalized_id_map[curie] = curie
                return curie

            # Molepro and probably others sometimes mistakenly use KEGG prefix
            #   instead of kEGG.COMPOUND... let biomapper choose between these
            if prefix.lower() == "kegg":
                prefix = ("kegg", "kegg.compound", "kegg.target")

            normalized_curie_dict, invalid_id_dict, unrecognized_vocabs = self.normalizer.get_curies(
                local_ids_dict={prefix: local_id}, stop_on_invalid_id=False, log_warnings=False, fuzzy_match_vocab=False
            )
            # Record curies it failed on
            self.unrecognized_vocabs |= unrecognized_vocabs
            if invalid_id_dict:
                for prefix, invalid_ids in invalid_id_dict.items():
                    self.prefixes_with_invalid_ids[prefix] += len(invalid_ids)

            # If it failed to normalize the curie, just return the original curie, unedited
            final_curie = list(normalized_curie_dict.keys())[0] if normalized_curie_dict else curie
            self.normalized_id_map[curie] = final_curie  # Cache our mapping
            return final_curie
        else:
            return curie
