import logging
from pathlib import Path

from kraken.harmonizers.base import BaseHarmonizer

# Trial records in the nodes file are identified by this prefix on their id (e.g. CLINICALTRIALS:NCT06519656).
TRIAL_PREFIX = "CLINICALTRIALS:"


class CTKGHarmonizer(BaseHarmonizer):
    """Clinical Trials KG (infores:multiomics-clinicaltrials): a TRAPI-style JSONL KG of
    intervention -> condition edges derived from clinical-trial data.

    Entity nodes are minimal ({id, name, category}); edges carry a TRAPI ``sources`` list (parsed
    by the base into primary/aggregator/supporting KS) plus clinical-trial context that
    flows into per-edge attributes.

    The nodes file also carries a second record shape -- one per clinical trial, with no category
    and a pile of clinical_trial_* metadata -- which we skip (see _stream_nodes)."""

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

    def _stream_nodes(self, input_path: Path | str):
        """Skip the clinical-trial records in the nodes file.

        The nodes file interleaves two record shapes: the entity nodes that edges actually refer to,
        and one record per trial (id CLINICALTRIALS:NCT..., no category, plus clinical_trial_* metadata
        and `interventions`/`conditions` lists). No edge has a trial as its subject or object -- trials
        appear only in the `has_supporting_studies` edge attribute, which keeps the NCT ids as-is -- and
        the interventions/conditions lists merely restate the intervention -> condition edges that are
        already present. So ingesting them would add only edgeless nodes, and we drop them here.

        Non-trial nodes are passed through untouched: one missing a category is a real problem with the
        source and should error out in _harmonize_node rather than be silently skipped."""
        skipped = 0
        for node in super()._stream_nodes(input_path):
            if node[self.id_prop].startswith(TRIAL_PREFIX):
                skipped += 1
                continue
            yield node
        if skipped:
            logging.info(
                f"Skipped {skipped} clinical-trial records in the {self.source_name} nodes file (no edge "
                f"refers to one; their NCT ids are retained in the 'has_supporting_studies' edge attribute)"
            )
