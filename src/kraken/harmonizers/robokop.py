from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import ROBOKOP_INFORES


class RobokopHarmonizer(BaseHarmonizer):
    source_name = "robokop"
    source_infores = ROBOKOP_INFORES
    is_aggregator = True

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = "equivalent_identifiers"
    synonyms_props = set()
    url_prop = "url"

    # Edge property config
    publications_info_prop = "sentences"
    primary_ks_exclusions = {"infores:ubergraph"}  # HUGE (60m edges) and we get it from Translator KG anyway
