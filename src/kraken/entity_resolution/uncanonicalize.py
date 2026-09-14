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

import re

from kraken.utils.constants import (
    EDGE_ATTRIBUTES,
    EDGE_OBJECT,
    EDGE_SUBJECT,
    ORIGINAL_OBJECT_ATTR,
    ORIGINAL_SUBJECT_ATTR,
)

# Sources whose stored endpoints are Babel-canonicalized (so match/edge remapping must
# recover the ORIGINAL endpoints). KG2 keeps its originals in ``kg2pre_ids`` (one merged
# edge can carry several); the other two store a single ``original_subject`` /
# ``original_object`` in their per-source attribute dict.
#
# The multiomics KGs (microbiome-kg, multiomics-kg) are deliberately NOT here, although they carry
# original_subject/original_object too. Theirs are not ids at all but the raw text of the paper tables they
# were built from ("schizophrenia", "Stool_sarcosine", "TAP1:Q03518:OID31289:v1") -- not one of their ~1.6M
# originals is a curie -- and their stored subject/object is their own curation of that text onto an id. So
# the stored endpoints are the authority. Listing them here remapped every edge onto a text label that never
# became a node: 2.1.1 kept 20 of multiomics-kg's 701,514 edges and 19 of microbiome-kg's 112,118.
CANONICALIZED_AGGREGATOR_SOURCES: frozenset[str] = frozenset({"kg2", "robokop", "translator-kg-open"})

# The subset whose recorded originals are also used as match-graph ids (aliases). kg2 is deliberately NOT here: all
# of kg2's equivalence is already in its equivalent_ids lists (only 168 of 312,596 sampled kg2 original ids are
# missing from them), and its originals are not reliably paired with the stored endpoints -- 12% of kg2 original
# pairs are SWAPPED relative to the stored edge (kg2 re-orients edges when it normalizes a relation), so pairing
# them positionally produced false aliases such as 80 DrugBank drugs "aliased" to one PathWhiz pathway. ROBOKOP's
# and Translator's originals are real ids nothing else supplies (HGVS -> CAID, Ensembl -> NCBIGene), 1:1 with the
# stored endpoint.
ALIAS_EVIDENCE_SOURCES: frozenset[str] = frozenset({"robokop", "translator-kg-open"})
# KG2's per-edge original ids: "orig_subject---relation---q---q---q---orig_object---src".
KG2_PRE_IDS_ATTR = "kg2pre_ids"
_KG2_ID_SEP = "---"
# (The other aggregators keep the original endpoints as plain attributes; those attribute names live
# in utils.constants because the harmonizers write them and entity resolution reads them.)


def kg2_pre_id_triples(edge: dict) -> list[tuple[str, str, str]]:
    """KG2's per-edge originals as (original_subject, original_relation, original_object).

    The relation matters for orientation: KG2 re-orients an edge when it normalizes an inverse relation
    (UMLS:RB "broader", NCIT:inverse_isa, HMDB:in_pathway, ...), so ~12% of original pairs run opposite to the
    stored edge -- and which relation it was decides that, with no exceptions in the data.
    """
    triples: list[tuple[str, str, str]] = []
    for attrs in (edge.get(EDGE_ATTRIBUTES) or {}).values():
        if not isinstance(attrs, dict):
            continue
        for raw in attrs.get(KG2_PRE_IDS_ATTR) or []:
            parts = raw.split(_KG2_ID_SEP)
            if len(parts) >= 6 and parts[0] not in ("", "None") and parts[5] not in ("", "None"):
                triples.append((parts[0], parts[1], parts[5]))
    return triples


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


def original_alias_pairs(edge: dict, source: str) -> list[tuple[str, str]]:
    """The ``(original_id, canonical_id)`` identity assertions this edge makes.

    When a canonicalizing aggregator stores an edge, it is recording "I resolved X to Y". That
    pairing is an equivalence claim of exactly the same kind as the ``equivalent_ids`` list on its
    nodes, so entity resolution consumes it the same way: as weighted, guardrail-checked match-graph
    evidence at the source's own equivalency weight -- never an unconditional merge. Doing so is
    what lets an edge keep its original endpoint as a real id instead of orphaning when that id
    exists nowhere in the node set (robokop's gtex edges are stored on ``CAID:`` nodes but originate
    at ``HGVS:`` ids; translator's originate at ``Ensembl:``/``MGI:``/``EMAPA:`` ids).

    Returns [] for a native source (its endpoints are already its own ids), for kg2 (see
    ALIAS_EVIDENCE_SOURCES), for an edge with no recoverable originals, and for any pair where the original
    and the canonical id are the same.
    """
    if source not in ALIAS_EVIDENCE_SOURCES:
        return []
    stored_subject, stored_object = edge.get(EDGE_SUBJECT), edge.get(EDGE_OBJECT)
    if not stored_subject or not stored_object:
        return []
    aliases: list[tuple[str, str]] = []
    for original_subject, original_object in original_endpoints(edge, source) or []:
        for original, canonical in ((original_subject, stored_subject), (original_object, stored_object)):
            if original and original != canonical:
                aliases.append((original, canonical))
    return aliases


# An original endpoint can name a COARSER entity than the node the aggregator canonicalized it to. The case that
# occurs: GWAS Catalog reports some associations against a bare dbSNP rsid -- a POSITION -- and ROBOKOP stores the
# edge on one ClinGen allele (CAID) at it. Taking that pairing as an alias would merge the position into that one
# allele (the same conflation ROBOKOP's own rsid equivalences caused; see harmonizers/robokop.py). An rsid WITH its
# allele ("DBSNP:rs142570322-T") names the allele, so it is not a mismatch and still aliases.
_POSITION_LEVEL_ORIGINAL = re.compile(r"DBSNP:rs\d+")
_ALLELE_LEVEL_PREFIX = "CAID:"


def is_coarser_than_canonical(original: str, canonical: str) -> bool:
    """True when ``original`` names a coarser entity than ``canonical`` (a bare rsid position vs. a CAID allele).

    Such a pair must not be used as equivalence evidence. The original is still a real id -- the edge should
    remap onto it, at the granularity the source actually asserted -- so callers keep seeding it.
    """
    return bool(_POSITION_LEVEL_ORIGINAL.fullmatch(original)) and canonical.startswith(_ALLELE_LEVEL_PREFIX)
