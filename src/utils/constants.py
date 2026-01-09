from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]

# Define node property names
ID = 'id'
NAME = 'name'
IRI = 'iri'
CATEGORIES = 'categories'
PROVIDED_BY = 'provided_by'
SYNONYMS = 'synonyms'
EQUIVALENT_IDS = 'equivalent_ids'
DESCRIPTION = 'description'
CHEMICAL_FORMULA = 'chemical_formula'
EXACT_MASS = 'exact_mass'
ATTRIBUTES = 'attributes'

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
QUALIFIED_DIRECTION = 'qualified_object_direction'
QUALIFIED_ASPECT = 'qualified_object_aspect'
CONTEXT_QUALIFIER = 'context_qualifier'

CORE_EDGE_PROPERTIES = {ID, SUBJECT, OBJECT, PREDICATE, PRIMARY_KS, AGGREGATOR_KS, SUPPORTING_SOURCES,
                        KNOWLEDGE_LEVEL, AGENT_TYPE, QUALIFIED_PREDICATE, QUALIFIED_DIRECTION, QUALIFIED_ASPECT}

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
ROBOKOP_INFORES: str = f"{INFORES_PREFIX}:robokop"
UMLS_INFORES: str = f'{INFORES_PREFIX}:umls'
LIPIDMAPS_CURIE: str = 'lipidmaps'
REFMET_CURIE: str = 'refmet'
CLINGEN_CURIE = f'{INFORES_PREFIX}:clingen'

KNOWN_INVALID = 'KNOWN_INVALID'


