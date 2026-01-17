from kraken.utils.constants import MICROBIOME_KG_INFORES
from kraken.harmonizers.base import BaseHarmonizer


class MicrobiomeKGHarmonizer(BaseHarmonizer):
    source_name = "microbiome-kg"
    source_infores = MICROBIOME_KG_INFORES

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = ""
    synonyms_props = set()
    url_prop = ""

    # Edge property config
    publications_prop = "publication"
    primary_ks_default_value = MICROBIOME_KG_INFORES
    supporting_sources_default_value = "infores:pubmed-central"
