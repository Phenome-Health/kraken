from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import TRANSLATOR_SOURCE_ID


class TranslatorKGOpenHarmonizer(BaseHarmonizer):
    source_infores = TRANSLATOR_SOURCE_ID
    is_aggregator = True

    # Node property config
    equivalent_ids_prop = "equivalent_identifiers"
