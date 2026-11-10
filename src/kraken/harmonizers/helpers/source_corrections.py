"""
Hand-checked corrections for single wrong statements inside a source's own records.

Each entry undoes one statement a source makes that contradicts the rest of that same record (its structure, its
other ids, its name), checked by hand. Entries are keyed by source and so touch only that source's records: the same
id from another source is left alone. A dropped id or synonym is not attached anywhere else, since the record never
said where it belongs. So an entry must not drop an id that has no record of its own: that id would stop being a
node, and every edge it is the original endpoint of would be dropped with it. The one exception is an id the source
lists on the wrong one of its own records: it is moved to the right record (MOVED_EQUIVALENT_IDS), so it keeps a
record and its edges follow it.

Edge corrections apply to sources harmonized from KGX files (``_harmonize_edges``); the single-file harmonizers write
their edges themselves. A source can fix its record in a later release, so the harmonizer logs every statement here
that it did not find; that is the signal to retire the entry.
"""

from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import NamedTuple


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
    # A screening-library macrocycle, "LSM-39319", carrying cotinine's synonym and description, which came with
    # PathWhiz's cotinine (moved below) (#29)
    ("infores:rtx-kg2", "CHEBI:127762"): NodeCorrection(
        synonyms=frozenset({"Cotinine"}),
        description_containing="the most important metabolite of nicotine",
    ),
}


class Move(NamedTuple):
    wrong_record: str
    right_record: str


# (source infores, id) -> the source's own records it moves between
MOVED_EQUIVALENT_IDS: dict[tuple[str, str], Move] = {
    # PathWhiz's cotinine, whose 8 nicotine pathway edges belong on cotinine, not on the macrocycle (#29)
    ("infores:rtx-kg2", "PathWhiz.Compound:824"): Move(wrong_record="CHEBI:127762", right_record="CHEBI:68641"),
}
# How the harmonizer records that a move found its right record. Missing, the moved id has no record at all.
MOVE_TARGET = "move_target"


def _moved_ids_by(record_of: Callable[[Move], str]) -> dict[tuple[str, str], frozenset[str]]:
    by_record: dict[tuple[str, str], set[str]] = defaultdict(set)
    for (source, curie), move in MOVED_EQUIVALENT_IDS.items():
        by_record[(source, record_of(move))].add(curie)
    return {record: frozenset(curies) for record, curies in by_record.items()}


# (source infores, record id) -> the ids moved off / onto that record
MOVED_OFF = _moved_ids_by(lambda move: move.wrong_record)
MOVED_ON = _moved_ids_by(lambda move: move.right_record)

# (source infores, subject, predicate, object, primary knowledge source), normalized, of edges to drop
EDGE_CORRECTIONS: frozenset[tuple[str, str, str, str, str]] = frozenset(
    {
        # Reactome maps its phenylacetate (phenylacetic acid, conjugated with glutamine in the mitochondria) to
        # phenyl acetate, the ester of phenol. Its own record keeps it a node (#29)
        ("infores:rtx-kg2", "REACT:R-ALL-177153", "biolink:same_as", "CHEBI:8082", "infores:reactome"),
        # HMDB's cotinine record (HMDB0001046, here under kg2's canonical id) lists the macrocycle's ChEBI id. A replay
        # of clustering shows this link alone pulls that ChEBI id into cotinine, away from its own InChIKey, once
        # PathWhiz's cotinine moves (#29)
        ("infores:rtx-kg2", "CHEBI:182609", "biolink:same_as", "CHEBI:127762", "infores:hmdb"),
    }
)

# (source infores, id, upstream source): drop the source's original statements in which that upstream source names
# the id, because it uses the id for something else. An edge left with no original statement is dropped. Only kg2
# keeps original statements, on each edge (kg2pre_ids); the id and the upstream source tag ("PathWhiz:") are spelled
# as kg2 writes them there, not normalized.
ORIGINAL_STATEMENT_CORRECTIONS: frozenset[tuple[str, str, str]] = frozenset(
    {
        # PathWhiz writes its cotinine as the macrocycle's ChEBI id too. Each of its 8 edges also names
        # PathWhiz.Compound:824, which the move puts on cotinine (#29)
        ("infores:rtx-kg2", "CHEBI:127762", "PathWhiz:"),
    }
)


def expected_corrections(source_infores: str) -> set[tuple[str, ...]]:
    """Every statement this source's corrections change, one tuple each, as the harmonizer records them applied."""
    expected: set[tuple[str, ...]] = {
        key for key in EDGE_CORRECTIONS | ORIGINAL_STATEMENT_CORRECTIONS if key[0] == source_infores
    }
    for (source, node_id), correction in NODE_CORRECTIONS.items():
        if source == source_infores:
            expected |= {(source, node_id, "equivalent_id", curie) for curie in correction.equivalent_ids}
            expected |= {(source, node_id, "synonym", synonym) for synonym in correction.synonyms}
            if correction.description_containing:
                expected.add((source, node_id, "description"))
    for (source, curie), move in MOVED_EQUIVALENT_IDS.items():
        if source == source_infores:
            expected.add((source, move.wrong_record, "equivalent_id", curie))
            expected.add((source, move.right_record, MOVE_TARGET, curie))
    return expected
