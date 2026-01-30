from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import KG2_INFORES


class KG2Harmonizer(BaseHarmonizer):
    source_name = "rtx-kg2"
    source_infores = KG2_INFORES
    is_aggregator = True

    # Node property config
    category_prop = "all_categories"
    equivalent_ids_prop = "equivalent_curies"
    synonyms_props = {"all_names"}
    url_prop = "iri"
    rename_node_attrs = {"category": "canonical_category"}

    # Edge property config
    ignore_edge_props = {"domain_range_exclusion"}
    rename_edge_attrs_or_quals = {
        "id": "kg2c_ids",
        "kg2_ids": "kg2pre_ids",
        "qualified_object_direction": "object_direction_qualifier",
        "qualified_object_aspect": "object_aspect_qualifier",
    }
    primary_ks_exclusions = {"infores:semmeddb"}
    # KG2 uses an invalid predicate for NCIT 'regimen_has_accepted_use_for_disease' edges - remap those
    predicate_overrides = {"biolink:drug_regulatory_status_world_wide": "biolink:treats_or_applied_or_studied_to_treat"}
