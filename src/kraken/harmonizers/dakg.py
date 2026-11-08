from typing import Any

from kraken.harmonizers.base import BaseHarmonizer

# What DAKG 0.5.7's own records still call the source, on every edge. The registry replaced it with
# infores:drugapprovals-kp, the source_id in build_config.
DEPRECATED_SELF_ID = "infores:multiomics-drugapprovals"


class DAKGHarmonizer(BaseHarmonizer):
    """Drug Approvals KG (infores:drugapprovals-kp): a TRAPI-style JSONL KG of
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

    def _parse_trapi_sources(self, edge: dict[str, Any]) -> tuple[str, list[str], list[str]]:
        # The raw `sources` list stays verbatim in the edge's attributes; only the flattened provenance is renamed.
        primary_ks, aggregator_kses, supporting_sources = super()._parse_trapi_sources(edge)
        current = {DEPRECATED_SELF_ID: self.source_infores}
        return (
            current.get(primary_ks, primary_ks),
            [current.get(aggregator, aggregator) for aggregator in aggregator_kses],
            supporting_sources,
        )
