from pathlib import Path

from kraken.schema import EdgeModel, NodeModel

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
UMLS_MTH_INFORES: str = f"{INFORES_PREFIX}:umls-metathesaurus"
REFMET_INFORES: str = f"{INFORES_PREFIX}:refmet"
CLINGEN_INFORES = f"{INFORES_PREFIX}:clingen"
HMDB_INFORES: str = f"{INFORES_PREFIX}:hmdb"


# Some sources do not (yet) have registered infores curies; use bare source IDs for those
NIH_CDE_SOURCE_ID: str = "nih-cde"
TRANSLATOR_SOURCE_ID: str = "translator-kg-open"
LIPIDMAPS_ID: str = "lipidmaps"
BIOLOGICAL_BMI_SOURCE_ID: str = "biological-bmi"  # multiomic BMI models, Watanabe et al. Nat Med 2023 (a paper, not a DB)
BIO_AGE_SOURCE_ID: str = "biological-age"  # multiomic biological-age models, Earls et al. J Gerontol A 2019 (a paper, not a DB)
PGS_CATALOG_SOURCE_ID: str = f"pgs-catalog"

# Primary knowledge sources whose edge publication lists are unreliable, so we drop publications from their
# edges during harmonization. (HMDB copies a disease's entire reference list onto every metabolite it links to
# that disease, so the PMIDs on an HMDB metabolite-disease edge describe the disease, not the specific
# assertion.) This is a temporary data-quality workaround; extend the set as more such sources are found.
UNRELIABLE_PUBLICATION_PRIMARY_KS: set[str] = {HMDB_INFORES}

KNOWN_INVALID = "KNOWN_INVALID"

NOT_PROVIDED = "not_provided"
MANUAL_AGENT = "manual_agent"
KNOWLEDGE_ASSERTION = "knowledge_assertion"
# Biolink KLAT values. For edges that report direct, dataset-specific statistical results (e.g. a feature's
# association with an outcome in a model's cohort), statistical_association pairs with data_analysis_pipeline.
# computational_model is for agents that generate broader conclusions/predictions (kept for such future edges).
STATISTICAL_ASSOCIATION = "statistical_association"
DATA_ANALYSIS_PIPELINE = "data_analysis_pipeline"
COMPUTATIONAL_MODEL = "computational_model"

NONE_STRINGS = {"none", "null", "-", "na", "n/a"}

QUALIFIED_PREDICATE = "qualified_predicate"
OBJ_DIRECTION_QUALIFIER = "object_direction_qualifier"
OBJ_ASPECT_QUALIFIER = "object_aspect_qualifier"

# ------------ For quick lookup during massive ETL! ------------- #

# Node property name constants
NODE_ID = NodeModel.id.name
NODE_NAME = NodeModel.name.name
NODE_URLS = NodeModel.urls.name
NODE_CATEGORIES = NodeModel.categories.name
NODE_PROVIDED_BY = NodeModel.provided_by.name
NODE_SYNONYMS = NodeModel.synonyms.name
NODE_EQUIVALENT_IDS = NodeModel.equivalent_ids.name
NODE_TAXA = NodeModel.taxa.name
NODE_DESCRIPTION = NodeModel.description.name
NODE_CHEMICAL_FORMULA = NodeModel.chemical_formula.name
NODE_EXACT_MASS = NodeModel.exact_mass.name
NODE_PUBLICATIONS = NodeModel.publications.name
NODE_ATTRIBUTES = NodeModel.attributes.name

# Edge property name constants
EDGE_SUBJECT = EdgeModel.subject.name
EDGE_OBJECT = EdgeModel.object.name
EDGE_PREDICATE = EdgeModel.predicate.name
EDGE_PRIMARY_KS = EdgeModel.primary_ks.name
EDGE_AGGREGATOR_KS = EdgeModel.aggregator_ks.name
EDGE_SUPPORTING_SOURCES = EdgeModel.supporting_sources.name
EDGE_KNOWLEDGE_LEVEL = EdgeModel.knowledge_level.name
EDGE_AGENT_TYPE = EdgeModel.agent_type.name
EDGE_QUALIFIERS = EdgeModel.qualifiers.name
EDGE_PUBLICATIONS = EdgeModel.publications.name
EDGE_PUBLICATIONS_INFO = EdgeModel.publications_info.name
EDGE_ATTRIBUTES = EdgeModel.attributes.name
