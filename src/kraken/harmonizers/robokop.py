from kraken.harmonizers.base import BaseHarmonizer


class RobokopHarmonizer(BaseHarmonizer):
    is_aggregator = True

    # Node property config
    category_prop = "category"
    equivalent_ids_prop = "equivalent_identifiers"
    synonyms_props = set()
    url_prop = "url"

    # Edge property config
    publications_info_prop = "sentences"
    # Not a source we ingest, so it stays a manual exclusion. HUGE (60m edges) and we get it from
    # Translator KG anyway -- which means Translator must actually be producing edges for this to
    # be a trade rather than a loss (it silently was not, in 2.1.1).
    source_exclusions = {"infores:ubergraph"}
