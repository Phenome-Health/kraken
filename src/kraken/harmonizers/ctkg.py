from kraken.harmonizers.base import BaseHarmonizer


class CTKGHarmonizer(BaseHarmonizer):
    """Clinical Trials KG (infores:multiomics-clinicaltrials): a TRAPI-style JSONL KG of
    intervention -> condition edges derived from clinical-trial data.

    Nodes are minimal ({id, name, category}); edges carry a TRAPI ``sources`` list (parsed
    by the base into primary/aggregator/supporting KS) plus clinical-trial context that
    flows into per-edge attributes."""

    # Nodes carry only id/name/category -- no xrefs, synonyms, or urls.
    equivalent_ids_prop = ""
    synonyms_props = set()
    url_prop = ""

    # Trust the TRAPI `sources` roles as given; don't force-record the KG itself as an aggregator.
    is_aggregator = False

    # Drop redundant / non-schema edge fields. The meaningful clinical-trial context
    # (max_research_phase, has_supporting_studies, elevate_to_prediction, tested_intervention,
    # intervention_boxed_warning) is retained automatically as per-edge attributes.
    ignore_edge_props = {"id", "category", "subject_name", "object_name"}
