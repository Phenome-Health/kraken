from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import MULTIOMICS_KG_INFORES


class MultiomicsKGHarmonizer(BaseHarmonizer):
    source_name = "multiomics-kg"
    source_infores = MULTIOMICS_KG_INFORES

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = ""
    synonyms_props = set()
    url_prop = ""

    # Edge property config
    publications_prop = "publication"
    primary_ks_default_value = MULTIOMICS_KG_INFORES
    supporting_sources_default_value = "infores:pubmed-central"
