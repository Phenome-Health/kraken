from kraken.harmonizers.base import BaseHarmonizer


class DAKGHarmonizer(BaseHarmonizer):
    """Drug Approvals KG (infores:multiomics-drugapprovals): a TRAPI-style JSONL KG of
    drug -> condition approval edges.

    Nodes are minimal ({id, name, category}); edges carry a TRAPI ``sources`` list (parsed
    by the base into primary/aggregator/supporting KS) plus approval context (N_cases,
    clinical_approval_status) that flows into per-edge attributes."""

    # Nodes carry only id/name/category -- no xrefs, synonyms, or urls.
    equivalent_ids_prop = ""
    synonyms_props = set()
    url_prop = ""

    # Trust the TRAPI `sources` roles as given; don't force-record the KG itself as an aggregator.
    is_aggregator = False

    # Drop redundant / non-schema edge fields; approval context (N_cases, clinical_approval_status)
    # is retained automatically as per-edge attributes.
    ignore_edge_props = {"id", "category", "subject_name", "object_name"}
