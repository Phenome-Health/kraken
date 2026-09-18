from kraken.config import get_source_id
from kraken.harmonizers.base import BaseHarmonizer


class MicrobiomeKGHarmonizer(BaseHarmonizer):
    # Node property config
    category_prop = "category"
    equivalent_ids_prop = ""
    synonyms_props = set()
    url_prop = ""

    # Edge property config
    publications_prop = "publication"
    supporting_sources_default_value = "infores:pubmed-central"
    primary_ks_default_value = get_source_id("microbiome-kg")
