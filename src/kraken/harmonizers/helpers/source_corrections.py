"""
Hand-checked corrections for single wrong statements inside a source's own records.

Each entry removes one statement a source makes that contradicts the rest of that same record (its structure, its
other ids, its name), checked by hand. Entries are keyed by (source, id) and so touch only that source's record:
the same id from another source is left alone. Nothing is redirected -- a dropped id or synonym is not attached
anywhere else, since the record never said where it belongs. So an entry must not drop an id that has no record of
its own: that id would stop being a node, and every edge it is the original endpoint of would be dropped with it.

Edge corrections apply to sources harmonized from KGX files (``_harmonize_edges``); the single-file harmonizers write
their edges themselves. A source can fix its record in a later release, so the harmonizer logs every statement here
that it did not find; that is the signal to retire the entry.
"""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class NodeCorrection:
    equivalent_ids: frozenset[str] = field(default_factory=frozenset)
    synonyms: frozenset[str] = field(default_factory=frozenset)
    # The record's description is dropped while it contains this; a repaired one stays (and the entry shows as stale)
    description_containing: str | None = None


NO_CORRECTION = NodeCorrection()

# (source infores, node id) -> what to drop from that source's record of the node, ids as the harmonizer normalizes them
NODE_CORRECTIONS: dict[tuple[str, str], NodeCorrection] = {
    # RefMet's econazole nitrate row lists N-succinyl-phenylalanine's ChEBI id; every other field is econazole (#29)
    ("infores:refmet", "RM:0048715"): NodeCorrection(equivalent_ids=frozenset({"CHEBI:232965"})),
    # Aceglutamide aluminium is a salt of acetyl-glutamine, not N-acetylasparagine (#29)
    ("infores:rtx-kg2", "CHEBI:139582"): NodeCorrection(synonyms=frozenset({"ACEGLUTAMIDE ALUMINIUM"})),
    # 4-iodoanisole, described under a glycopeptide's name ("...diacetylchitobiosyl-L-asparagine, isoform B") (#29)
    ("infores:rtx-kg2", "CHEBI:174192"): NodeCorrection(
        description_containing="diacetylchitobiosyl-L-asparagine, isoform B"
    ),
}

# (source infores, subject, predicate, object, primary knowledge source), normalized, of edges to drop
EDGE_CORRECTIONS: frozenset[tuple[str, str, str, str, str]] = frozenset(
    {
        # Reactome maps its phenylacetate (phenylacetic acid, conjugated with glutamine in the mitochondria) to
        # phenyl acetate, the ester of phenol. Its own record keeps it a node (#29)
        ("infores:rtx-kg2", "REACT:R-ALL-177153", "biolink:same_as", "CHEBI:8082", "infores:reactome"),
    }
)


def expected_corrections(source_infores: str) -> set[tuple[str, ...]]:
    """Every statement this source's corrections remove, one tuple each, as the harmonizer records them applied."""
    expected: set[tuple[str, ...]] = {key for key in EDGE_CORRECTIONS if key[0] == source_infores}
    for (source, node_id), correction in NODE_CORRECTIONS.items():
        if source == source_infores:
            expected |= {(source, node_id, "equivalent_id", curie) for curie in correction.equivalent_ids}
            expected |= {(source, node_id, "synonym", synonym) for synonym in correction.synonyms}
            if correction.description_containing:
                expected.add((source, node_id, "description"))
    return expected
