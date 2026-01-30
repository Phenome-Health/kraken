from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[3]


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
