from pathlib import Path

# Define node property names
ID = 'id'
NAME = 'name'
URLS = 'urls'
CATEGORIES = 'categories'
PROVIDED_BY = 'provided_by'
SYNONYMS = 'synonyms'
EQUIVALENT_IDS = 'equivalent_ids'
DESCRIPTION = 'description'
CHEMICAL_FORMULA = 'chemical_formula'
EXACT_MASS = 'exact_mass'
ATTRIBUTES = 'attributes'

CORE_NODE_PROPERTIES = {ID, NAME, URLS, CATEGORIES, PROVIDED_BY, SYNONYMS, EQUIVALENT_IDS,
                        DESCRIPTION, CHEMICAL_FORMULA, EXACT_MASS}

# Define edge property names
SUBJECT = 'subject'
OBJECT = 'object'
PREDICATE = 'predicate'
PRIMARY_KS = 'primary_knowledge_source'
AGGREGATOR_KS = 'aggregator_knowledge_source'
SUPPORTING_SOURCES = 'supporting_data_sources'
KNOWLEDGE_LEVEL = 'knowledge_level'
AGENT_TYPE = 'agent_type'
QUALIFIED_PREDICATE = 'qualified_predicate'
OBJ_DIRECTION_QUALIFIER = 'object_direction_qualifier'
OBJ_ASPECT_QUALIFIER = 'object_aspect_qualifier'
CONTEXT_QUALIFIER = 'context_qualifier'

CORE_EDGE_PROPERTIES = {ID, SUBJECT, OBJECT, PREDICATE, PRIMARY_KS, AGGREGATOR_KS, SUPPORTING_SOURCES,
                        KNOWLEDGE_LEVEL, AGENT_TYPE, QUALIFIED_PREDICATE, OBJ_DIRECTION_QUALIFIER, OBJ_ASPECT_QUALIFIER}

PUBLICATIONS = 'publications'
PUBLICATIONS_INFO = 'publications_info'


ROOT_CATEGORY = 'biolink:NamedThing'
ROOT_PREDICATE = 'biolink:related_to'

UNKNOWN_KNOWLEDGE_LEVEL = 'not_provided'
UNKNOWN_AGENT_TYPE = 'not_provided'

BIOLINK_PREFIX = 'biolink'
INFORES_PREFIX = 'infores'

SPOKE_INFORES: str = f'{INFORES_PREFIX}:spoke'
KG2_INFORES: str = f'{INFORES_PREFIX}:rtx-kg2'
ROBOKOP_INFORES: str = f"{INFORES_PREFIX}:robokop-kg"
MOLEPRO_INFORES: str = f"{INFORES_PREFIX}:molepro"
MICROBIOME_KG_INFORES: str = f"{INFORES_PREFIX}:multiomics-microbiome"
MULTIOMICS_KG_INFORES: str = f"{INFORES_PREFIX}:multiomics-multiomics"
UMLS_INFORES: str = f'{INFORES_PREFIX}:umls'
REFMET_INFORES: str = f'{INFORES_PREFIX}:refmet'
CLINGEN_INFORES = f'{INFORES_PREFIX}:clingen'

LIPIDMAPS_ID: str = 'lipidmaps'

KNOWN_INVALID = 'KNOWN_INVALID'

NOT_PROVIDED = 'not_provided'

NONE_STRINGS = {"none", "null", "-", "na", "n/a"}


