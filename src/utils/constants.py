

# Define core node property names
ID = 'id'
NAME = 'name'
IRI = 'iri'
CATEGORIES = 'categories'
PROVIDED_BY = 'provided_by'
SYNONYMS = 'synonyms'
EQUIVALENT_IDS = 'equivalent_ids'

CORE_NODE_PROPERTIES = {ID, NAME, IRI, CATEGORIES, PROVIDED_BY, SYNONYMS, EQUIVALENT_IDS}

CHEMICAL_FORMULA = 'chemical_formula'
EXACT_MASS = 'exact_mass'



# Define core edge property names
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

CORE_EDGE_PROPERTIES = {ID, SUBJECT, OBJECT, PREDICATE, PRIMARY_KS, AGGREGATOR_KS, SUPPORTING_SOURCES,
                        KNOWLEDGE_LEVEL, AGENT_TYPE, QUALIFIED_PREDICATE, QUALIFIED_DIRECTION, QUALIFIED_ASPECT}


ROOT_CATEGORY = 'biolink:NamedThing'
ROOT_PREDICATE = 'biolink:related_to'

UNKNOWN_KNOWLEDGE_LEVEL = 'not_provided'
UNKNOWN_AGENT_TYPE = 'not_provided'

SPOKE_INFORES: str = 'infores:spoke'
KG2_INFORES: str = 'infores:kg2'
UMLS_INFORES: str = 'infores:umls'
LIPIDMAPS_CURIE: str = 'lipidmaps'
REFMET_CURIE: str = 'refmet'

BIOLINK_PREFIX = 'biolink'
INFORES_PREFIX = 'infores'

KNOWN_INVALID = 'KNOWN_INVALID'
