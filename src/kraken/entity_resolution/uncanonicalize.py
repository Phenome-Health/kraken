"""Un-canonicalization of Babel-canonicalized aggregator edges.

A canonicalizing aggregator (kg2, robokop, ...) stores each edge on the Babel
canonical endpoints and keeps the ORIGINAL ids in an edge attribute. We must map
edges by their *original* endpoints -- both when building match-graph evidence and
when remapping the final edges -- because our clustering deliberately diverges from
Babel, so a canonical endpoint would attach an edge to the wrong node. One merged
edge can carry several original pairs (kg2), which become several edges.

Shared by ``build`` (match evidence) and ``integrate`` (edge remapping) so there is
exactly one un-canonicalization implementation.
"""

from __future__ import annotations

from kraken.utils.constants import EDGE_ATTRIBUTES, EDGE_OBJECT, EDGE_SUBJECT

# Sources whose stored endpoints are Babel-canonicalized (so match/edge remapping must
# recover the ORIGINAL endpoints). KG2 keeps its originals in ``kg2pre_ids`` (one merged
# edge can carry several); the other four store a single ``original_subject`` /
# ``original_object`` in their per-source attribute dict.
CANONICALIZED_AGGREGATOR_SOURCES: frozenset[str] = frozenset(
    {"kg2", "robokop", "translator-kg-open", "microbiome-kg", "multiomics-kg"}
)
# KG2's per-edge original ids: "orig_subject---relation---q---q---q---orig_object---src".
KG2_PRE_IDS_ATTR = "kg2pre_ids"
_KG2_ID_SEP = "---"
# The other aggregators keep the original endpoints as plain attributes.
ORIGINAL_SUBJECT_ATTR = "original_subject"
ORIGINAL_OBJECT_ATTR = "original_object"


def _kg2_pre_id_pairs(edge: dict) -> list[tuple[str, str]]:
    # kg2pre_ids lives in KG2's per-source attribute dict; find it regardless of the infores key it's under
    # (that key is KG2's build_config source_id, the single source of truth -- not restated here).
    pairs: list[tuple[str, str]] = []
    for attrs in (edge.get(EDGE_ATTRIBUTES) or {}).values():
        if not isinstance(attrs, dict):
            continue
        for raw in attrs.get(KG2_PRE_IDS_ATTR) or []:
            parts = raw.split(_KG2_ID_SEP)
            if len(parts) >= 6 and parts[0] not in ("", "None") and parts[5] not in ("", "None"):
                pairs.append((parts[0], parts[5]))
    return pairs


def _attribute_original_pair(edge: dict) -> tuple[str, str] | None:
    """The single (original_subject, original_object) an aggregator stores in its
    per-source attribute dict (found regardless of the infores key it's under)."""
    for attrs in (edge.get(EDGE_ATTRIBUTES) or {}).values():
        if isinstance(attrs, dict):
            subj, obj = attrs.get(ORIGINAL_SUBJECT_ATTR), attrs.get(ORIGINAL_OBJECT_ATTR)
            if subj and obj:
                return (subj, obj)
    return None


def original_endpoints(edge: dict, source: str) -> list[tuple[str, str]] | None:
    """Original (pre-canonicalization) subject/object pairs for an edge.

    Native (non-canonicalized) sources use their own subject/object. Canonicalizing
    aggregators recover the originals: KG2 from ``kg2pre_ids`` (possibly several pairs
    per merged edge); the rest from a single ``original_subject`` / ``original_object``
    attribute. Returns ``None`` only if a canonicalized aggregator edge carries no
    recoverable originals, so the caller decides the fallback (skip a match edge; use
    the stored canonical endpoints when remapping a real edge).
    """
    if source not in CANONICALIZED_AGGREGATOR_SOURCES:
        return [(edge.get(EDGE_SUBJECT, ""), edge.get(EDGE_OBJECT, ""))]
    if source == "kg2":
        pairs = _kg2_pre_id_pairs(edge)
        if pairs:
            return pairs
    single = _attribute_original_pair(edge)  # robokop/translator/mokg/mbkg (and kg2 w/o pre_ids)
    return [single] if single else None
