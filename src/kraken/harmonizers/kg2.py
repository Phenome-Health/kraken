from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import KG2_INFORES


class KG2Harmonizer(BaseHarmonizer):
    source_name = "rtx-kg2"
    source_infores = KG2_INFORES

    # Node property config
    category_prop = "all_categories"
    equivalent_ids_prop = "equivalent_curies"
    synonyms_props = {"all_names"}
    url_prop = "iri"
    rename_node_attrs = {"category": "canonical_category"}

    # Edge property config
    object_direction_qualifier_prop = "qualified_object_direction"
    object_aspect_qualifier_prop = "qualified_object_aspect"
    ignore_edge_props = {"domain_range_exclusion"}
    rename_edge_attrs = {"id": "kg2c_ids", "kg2_ids": "kg2pre_ids"}
    # KG2 uses an invalid predicate for NCIT 'regimen_has_accepted_use_for_disease' edges - remap those
    predicate_overrides = {"biolink:drug_regulatory_status_world_wide": "biolink:treats_or_applied_or_studied_to_treat"}
