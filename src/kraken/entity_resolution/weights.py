"""Per-source / per-predicate evidence weights for the match graph.

Design (see ``docs/entity_resolution_plan.md`` §1):

* Weights are **per source** and **accumulate across agreeing sources**.
* Clustering is label propagation (no resolution parameter), so the gate is
  **tau**: a CURIE pair whose accumulated weight reaches tau gets an edge, and
  label propagation then merges what is connected (subject to the guardrails). So
  a source whose weight alone reaches tau can merge a pair on its own; a source
  below tau can only *corroborate* (its weight sums with other sub-tau evidence
  to cross tau). Only ``close_match`` sits below tau (a weak "roughly the same"
  assertion); every equivalence source — including Babel's cliques and the aggregators'
  Babel-derived lists — sits **at or above** tau so it can merge on its own. The aggregators
  are kept safe not by a sub-tau weight but by (a) sharing one source group
  (max, not sum) so echoing Babel counts once, (b) the guardrail edge-prune,
  which drops cross-family conflation edges before clustering regardless of weight,
  and (c) being **prefix-capped**: ids of a prefix a list holds in bulk don't merge on
  their own, because that bulk is where aggregator lists go wrong (see ``max_ids_per_prefix``).
* Correlated sources are de-correlated: KG2 / ROBOKOP / Translator all derive
  equivalence from the SRI Node Normalizer, i.e. from Babel, so neither they nor Babel itself are
  independent evidence. Evidence within a source group is combined by
  **max**, not sum, so their shared ancestry cannot triple-count and re-import
  Babel's clustering. Independent sources accumulate by sum.

ALL NUMBERS HERE ARE UNTUNED PLACEHOLDERS. They exist so the pipeline runs
end-to-end; the eval harness (``eval/scorer.py``) is what tunes them. Override
via ``config/entity_resolution/weights.yaml``.
"""

from __future__ import annotations

from pathlib import Path

import yaml
from pydantic import BaseModel, Field

from kraken.utils.constants import CLOSE_MATCH_PREDICATE, PROJECT_ROOT, SAME_AS_PREDICATE

# Optional tuning file; if absent, the defaults below apply.
DEFAULT_WEIGHTS_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "weights.yaml"

# Match predicates. exact_match / same_as are full-strength equivalence;
# close_match is weak; broad_match / narrow_match are hierarchical and MUST be
# excluded (including them guarantees parent/child collapse). Predicate values
# are always biolink-prefixed after harmonization, so only prefixed forms appear.
EXACT_MATCH_PREDICATES: frozenset[str] = frozenset({"biolink:exact_match", SAME_AS_PREDICATE})
CLOSE_MATCH_PREDICATES: frozenset[str] = frozenset({CLOSE_MATCH_PREDICATE})
EXCLUDED_MATCH_PREDICATES: frozenset[str] = frozenset({"biolink:broad_match", "biolink:narrow_match"})

NAME_SIMILARITY_GROUP = "name_similarity"


class ERWeights(BaseModel):
    """Weights + thresholds for match-graph construction and clustering."""

    # Edge-existence threshold. A CURIE pair whose accumulated weight reaches tau
    # gets an edge (and label propagation will then merge it, guardrails allowing);
    # below tau it is dropped. This is THE merge gate under label propagation.
    tau: float = 0.3

    # Per-source weight for an equivalency-list assertion (one clique edge).
    # Sources absent from this map use ``default_equivalency_weight``.
    equivalency_weights: dict[str, float] = Field(
        default_factory=lambda: {
            # THE equivalence backbone: Babel's cliques and gene/protein conflations, read from its release
            # files (harmonizers/babel.py). Above tau so a clique merges on its own; cross-family conflations
            # are pruned pairwise before clustering, so this only governs SAME-family merges.
            "babel": 0.5,
            # Curated, structurally tight -> reach tau, so they merge on their own.
            "ncbigene": 1.0,
            "refmet": 1.0,
            "lipidmaps": 1.0,
            "umls": 0.8,
            "loinc": 0.6,
            "cdes": 0.6,
            # Aggregators' baked-in equivalent_ids lists. At or above tau, so a list can merge on
            # its own -- except ids of a prefix it holds in bulk: these sources are PREFIX-CAPPED
            # (see max_ids_per_prefix), which is what actually keeps them safe. Their value is
            # recovering the mappings Babel doesn't know, and that gap is large: of kg2's list pairs
            # 64% were unknown to the Node Normalizer (RXCUI/UMLS/RXNORM/CHV), of robokop's 78% (almost all DBSNP).
            # Sharing the "babel_derived" group (max-not-sum) still stops the same Babel
            # assertion counting twice across aggregators.
            #
            # (These were a flat 0.15 -- below tau, and with nothing else sub-tau in the graph to
            # sum with, that made them inert: no list ever created an edge. It is why every
            # robokop DBSNP id shipped in 2.1.1 as a nameless singleton beside its CAID node.)
            "kg2": 0.5,
            "robokop": 0.5,
            "translator-kg-open": 0.5,
        }
    )
    default_equivalency_weight: float = 0.4

    # PREFIX-CAPPED equivalency. An aggregator's equivalent_ids list goes wrong in one characteristic way:
    # ONE prefix contributing a pile of ids that aren't the entity -- kg2 lists 1,329 Reactome reactions on TP53,
    # 132 RxCUI branded products on metformin -- while the prefixes contributing a few ids (TP53's NCIT gene
    # concept, a disease's CHV and OMIM:MTHU terms, a drug's ATC code) are good. Judging the whole list by its
    # size threw out the good ids with the bulk: every one of those lost ids came only from kg2 and ended up a
    # singleton.
    #
    # So within a list, ids whose prefix appears more than ``max_ids_per_prefix[source]`` times are BULK: each
    # gets a single ``bulk_prefix_weight`` edge to the node that listed it (below tau -- it records that they are
    # related, and can add to other evidence, but never merges on its own; a star, so 1,329 bulk ids cost 1,329
    # edges rather than ~880k). Every other id keeps the source's full weight, however long the list is.
    #
    # Only the sources listed are capped. Everything else is trusted as a whole list -- above all "babel", whose
    # cliques are clean and legitimately hold many ids of one prefix (a gene's protein isoforms). The cap is a
    # starting value: it clears the bulk seen so far (hundreds per prefix), and nothing measurable separates good
    # from bad at 11-50 ids of one prefix, so tune it against spot checks.
    max_ids_per_prefix: dict[str, int] = Field(
        default_factory=lambda: {
            "kg2": 10,
            "robokop": 10,
            "translator-kg-open": 10,
        }
    )
    bulk_prefix_weight: float = 0.05
    #
    # An ALIAS -- the original -> canonical pairing a canonicalized aggregator records on an edge
    # (see uncanonicalize.original_alias_pairs) -- is the same kind of assertion as a list entry
    # ("this source says these are one thing"), so it takes the source's full equivalency weight rather than
    # a knob of its own, and the two can never drift apart.

    # Per-predicate weights for source match-predicate edges. close_match is a weak
    # "roughly the same" assertion (not true equivalence), so it's kept very low —
    # well below tau, so it only ever corroborates, never merges on its own (and the
    # subclass co-occurrence penalty can drive it lower still).
    #
    # Corroboration here means PARALLEL assertions: several close_match edges on the same pair, each
    # from a different primary knowledge source, are independent claims and should sum -- three at
    # 0.1 reach tau. So match-predicate evidence is grouped per (source, primary knowledge source),
    # not per source (see predicate_group). Grouping by source alone put every kg2 assertion in
    # babel_derived, where they combine by MAX and three parallel edges came to 0.1, not 0.3.
    #
    # Today this rarely fires: on the ORIGINAL endpoints ER matches on, no kg2 pair is asserted by
    # two different KSes. What looks like agreement on kg2's stored endpoints is different source
    # ids (ATC, CHV, MeSH, ...) that Babel folded onto one node -- counting that would re-import
    # Babel's clustering, so un-canonicalizing is correct not to see it. The handling is here for
    # sources that do assert parallel matches on the same ids.
    exact_match_weight: float = 1.0
    close_match_weight: float = 0.1

    # Name-similarity edges: at/above tau, so two nodes linked ONLY by a shared
    # (normalized, primary) name still merge. This is safe because the guardrail
    # edge-prune at formation already removes name edges between incompatible
    # nodes (different branch / taxon / enforced structural id), so a name edge
    # that survives is branch-, taxon-, and structural-id-compatible — strong
    # evidence, not weak. The residual risk is two DISTINCT entities in the SAME
    # branch/taxon that happen to share a normalized name and carry no enforced
    # id (e.g. two CHEBI with identical labels); the eval measures that cost.
    # UNTUNED placeholder like the rest — the eval sets the final value.
    name_similarity_weight: float = 0.7

    # A close_match edge between two nodes that ALSO have subclass_of/superclass_of
    # edges between them is likely a mislabeled hierarchical relation, not
    # equivalence. Each such hierarchical edge multiplies the close_match weight by
    # this decay (0.5 -> one halves it, two quarters it, ...), so more hierarchical
    # evidence -> weaker close_match. 1.0 disables the penalty. UNTUNED.
    subclass_penalty_decay: float = 0.5

    # Source groups: sources listed together contribute by MAX, not sum.
    # Each source maps to a group id; sources not listed are their own group.
    source_groups: dict[str, list[str]] = Field(
        default_factory=lambda: {
            # All Babel-derived equivalence shares one group so the same Babel assertion
            # (Babel's own clique + the aggregators' baked-in lists) counts once (max), not summed.
            "babel_derived": ["babel", "kg2", "robokop", "translator-kg-open"],
        }
    )

    # Equivalency clique handling. Sets up to this size become full cliques; larger
    # sets become a star from the list's head. Now that equivalence comes
    # from Babel's CLEAN cliques (not aggregators' junk lists), this is a
    # pure SCALE valve (avoid N^2 edges on a pathological clique), not a distrust
    # mechanism — so keep it above real clique sizes (gene/protein ~30-40, disease
    # ~15-20) so legitimate cliques stay full (robust), and only true outliers star.
    clique_cap: int = 100

    # Name-similarity blocking: skip groups larger than this (likely a generic
    # token), and drop names shorter than this or purely numeric.
    name_group_cap: int = 40
    min_name_length: int = 3

    def model_post_init(self, _context: object) -> None:
        if self.bulk_prefix_weight >= self.tau:
            raise ValueError(
                f"bulk_prefix_weight ({self.bulk_prefix_weight}) must stay below tau ({self.tau}); "
                f"otherwise a list's bulk-prefix ids would merge on their own"
            )
        # Precompute source -> source-group-id.
        self._source_to_group: dict[str, str] = {}
        for group_id, members in self.source_groups.items():
            for member in members:
                self._source_to_group[member] = group_id

    # ---- lookups ----

    @property
    def aggregator_list_sources(self) -> set[str]:
        """Sources whose equivalency lists are aggregator lists: prefix-capped, and emitted as a STAR from the
        listing node rather than a clique (see match_graph.clique_evidence)."""
        return set(self.max_ids_per_prefix)

    def source_group(self, source: str) -> str:
        """Group id used for source-group de-correlation; independent sources get their own."""
        return self._source_to_group.get(source, f"src:{source}")

    def predicate_group(self, source: str, primary_ks: str | None) -> str:
        """Group id for a match-predicate edge: its source's group, split by primary knowledge source.

        Parallel close_match edges on one pair from DIFFERENT primary KSes are independent claims, so
        they land in different groups and sum (see close_match_weight). An edge naming no primary KS
        falls back to the plain source group, so unattributed claims from one source still combine by
        max rather than counting as independent.
        """
        group = self.source_group(source)
        return f"{group}|ks:{primary_ks}" if primary_ks else group

    def equivalency_weight(self, source: str) -> float:
        """Weight for one clique edge from this source's equivalency lists (bulk-prefix ids excepted -- see
        ``max_ids_per_prefix``)."""
        return self.equivalency_weights.get(source, self.default_equivalency_weight)

    def alias_weight(self, source: str) -> float:
        """Weight for one original -> canonical alias: the source's full equivalency weight."""
        return self.equivalency_weight(source)

    def predicate_weight(self, predicate: str) -> float | None:
        """Weight for a match-predicate edge, or ``None`` if the predicate must
        not contribute (hierarchical, or not a match predicate)."""
        if predicate in EXCLUDED_MATCH_PREDICATES:
            return None
        if predicate in EXACT_MATCH_PREDICATES:
            return self.exact_match_weight
        if predicate in CLOSE_MATCH_PREDICATES:
            return self.close_match_weight
        return None

    @classmethod
    def load(cls, path: str | Path | None = DEFAULT_WEIGHTS_PATH) -> ERWeights:
        """Load from YAML, falling back to defaults if the file is absent."""
        if path is None or not Path(path).exists():
            return cls()
        data = yaml.safe_load(Path(path).read_text()) or {}
        return cls(**data)
