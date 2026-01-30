from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[3]

# Define node property names
ID = "id"
NAME = "name"
URLS = "urls"
CATEGORIES = "categories"
PROVIDED_BY = "provided_by"
SYNONYMS = "synonyms"
EQUIVALENT_IDS = "equivalent_ids"
DESCRIPTION = "description"
CHEMICAL_FORMULA = "chemical_formula"
EXACT_MASS = "exact_mass"
ATTRIBUTES = "attributes"

CORE_NODE_PROPERTIES = {
    ID: str,
    NAME: str,
    URLS: list,
    CATEGORIES: list,
    PROVIDED_BY: list,
    SYNONYMS: list,
    EQUIVALENT_IDS: list,
    DESCRIPTION: str,
    CHEMICAL_FORMULA: str,
    EXACT_MASS: float,
}

# Define edge property names
SUBJECT = "subject"
OBJECT = "object"
PREDICATE = "predicate"
PRIMARY_KS = "primary_knowledge_source"
AGGREGATOR_KS = "aggregator_knowledge_source"
SUPPORTING_SOURCES = "supporting_data_sources"
KNOWLEDGE_LEVEL = "knowledge_level"
AGENT_TYPE = "agent_type"
QUALIFIERS = "qualifiers"
PUBLICATIONS = "publications"
PUBLICATIONS_INFO = "publications_info"

REQUIRED = "required"
TYPE = "type"
IN_KEY = "in-key"

# TODO: convert this to a little class and instantiate
EDGE_PROPERTIES = {
    SUBJECT: {REQUIRED: True, TYPE: str, IN_KEY: True},
    OBJECT: {REQUIRED: True, TYPE: str, IN_KEY: True},
    PREDICATE: {REQUIRED: True, TYPE: str, IN_KEY: True},
    PRIMARY_KS: {REQUIRED: True, TYPE: str, IN_KEY: True},
    SUPPORTING_SOURCES: {REQUIRED: False, TYPE: list, IN_KEY: True},
    QUALIFIERS: {REQUIRED: False, TYPE: dict, IN_KEY: True},
    KNOWLEDGE_LEVEL: {REQUIRED: True, TYPE: str, IN_KEY: False},
    AGENT_TYPE: {REQUIRED: True, TYPE: str, IN_KEY: False},
    AGGREGATOR_KS: {REQUIRED: False, TYPE: list, IN_KEY: True},
    PUBLICATIONS: {REQUIRED: False, TYPE: list, IN_KEY: True},
    PUBLICATIONS_INFO: {REQUIRED: False, TYPE: dict, IN_KEY: True},
    ATTRIBUTES: {REQUIRED: False, TYPE: dict, IN_KEY: True},
}
EDGE_KEY_PROPERTIES = {slot_name for slot_name, info in EDGE_PROPERTIES.items() if info[IN_KEY]}
REQUIRED_EDGE_PROPERTIES = {slot_name for slot_name, info in EDGE_PROPERTIES.items() if info[REQUIRED]}

ROOT_CATEGORY = "biolink:NamedThing"
ROOT_PREDICATE = "biolink:related_to"

BIOLINK_PREFIX = "biolink"
INFORES_PREFIX = "infores"


SPOKE_INFORES: str = f"{INFORES_PREFIX}:spoke"
KG2_INFORES: str = f"{INFORES_PREFIX}:rtx-kg2"
ROBOKOP_INFORES: str = f"{INFORES_PREFIX}:robokop-kg"
MOLEPRO_INFORES: str = f"{INFORES_PREFIX}:molepro"
MICROBIOME_KG_INFORES: str = f"{INFORES_PREFIX}:multiomics-microbiome"
MULTIOMICS_KG_INFORES: str = f"{INFORES_PREFIX}:multiomics-multiomics"
UMLS_INFORES: str = f"{INFORES_PREFIX}:umls"
REFMET_INFORES: str = f"{INFORES_PREFIX}:refmet"
CLINGEN_INFORES = f"{INFORES_PREFIX}:clingen"

LIPIDMAPS_ID: str = "lipidmaps"

KNOWN_INVALID = "KNOWN_INVALID"

NOT_PROVIDED = "not_provided"

NONE_STRINGS = {"none", "null", "-", "na", "n/a"}


QUALIFIED_PREDICATE = "qualified_predicate"
OBJ_DIRECTION_QUALIFIER = "object_direction_qualifier"
OBJ_ASPECT_QUALIFIER = "object_aspect_qualifier"
CONTEXT_QUALIFIER = "context_qualifier"
