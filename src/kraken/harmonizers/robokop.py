from .base import BaseHarmonizer
from kraken.utils.constants import ROBOKOP_INFORES


class RobokopHarmonizer(BaseHarmonizer):
    source_name = "robokop"
    source_infores = ROBOKOP_INFORES

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = "equivalent_identifiers"
    synonyms_props = set()  # TODO: check on this? are there any?
    url_prop = ""  # TODO: check on this? do they give iri or anything?

    # Edge property config
    publications_info_prop = "sentences"

