from dataclasses import dataclass


@dataclass(frozen=True)
class PropertyDef:
    """Metadata about a single property (not an instance of the property itself)."""

    name: str
    dtype: type | tuple[type, ...]
    required: bool = False
    inner_type: type | None = None
    in_key: bool = False


class NodeModel:
    """Defines the schema/properties for nodes (not a node instance)."""

    id = PropertyDef("id", str, required=True)
    name = PropertyDef("name", str)
    urls = PropertyDef("urls", list, inner_type=str)
    categories = PropertyDef("categories", list, required=True, inner_type=str)
    provided_by = PropertyDef("provided_by", list, required=True, inner_type=str)
    synonyms = PropertyDef("synonyms", list, inner_type=str)
    equivalent_ids = PropertyDef("equivalent_ids", list, inner_type=str)
    description = PropertyDef("description", str)
    chemical_formula = PropertyDef("chemical_formula", str)
    exact_mass = PropertyDef("exact_mass", float)
    publications = PropertyDef("publications", list, inner_type=str)
    attributes = PropertyDef("attributes", dict)

    @classmethod
    def all_properties(cls) -> dict[str, PropertyDef]:
        return {k: v for k, v in vars(cls).items() if isinstance(v, PropertyDef)}

    @classmethod
    def required_properties(cls) -> set[str]:
        return {p.name for p in cls.all_properties().values() if p.required}


class EdgeModel:
    """Defines the schema/properties for edges (not an edge instance)."""

    subject = PropertyDef("subject", str, required=True, in_key=True)
    object = PropertyDef("object", str, required=True, in_key=True)
    predicate = PropertyDef("predicate", str, required=True, in_key=True)
    qualifiers = PropertyDef("qualifiers", dict, in_key=True)
    primary_ks = PropertyDef("primary_knowledge_source", str, required=True, in_key=True)
    aggregator_ks = PropertyDef("aggregator_knowledge_source", list, in_key=True, inner_type=str)
    supporting_sources = PropertyDef("supporting_data_sources", list, in_key=True, inner_type=str)
    knowledge_level = PropertyDef("knowledge_level", str, required=True)
    agent_type = PropertyDef("agent_type", str, required=True)
    publications = PropertyDef("publications", list, inner_type=str)
    publications_info = PropertyDef("publications_info", (dict, str))
    attributes = PropertyDef("attributes", dict)

    @classmethod
    def all_properties(cls) -> dict[str, PropertyDef]:
        return {k: v for k, v in vars(cls).items() if isinstance(v, PropertyDef)}

    @classmethod
    def required_properties(cls) -> set[str]:
        return {p.name for p in cls.all_properties().values() if p.required}

    @classmethod
    def key_properties(cls) -> set[str]:
        return {p.name for p in cls.all_properties().values() if p.in_key}
