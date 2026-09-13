"""Per-source / per-predicate evidence weights for the match graph.

Design (see ``docs/entity_resolution_plan.md`` §1):

* Weights are **per source** and **accumulate across agreeing sources**.
* Clustering is label propagation (no resolution parameter), so the gate is
  **tau**: a CURIE pair whose accumulated weight reaches tau gets an edge, and
  label propagation then merges what is connected (subject to the guardrails). So
  a source whose weight alone reaches tau can merge a pair on its own; a source
  below tau can only *corroborate* (its weight sums with other sub-tau evidence
  to cross tau). Only ``close_match`` sits below tau (a weak "roughly the same"
  assertion); every equivalence source — including the aggregators' SRI/Babel
  cliques — sits **at or above** tau so it can merge on its own. The aggregators
  are kept safe not by a sub-tau weight but by (a) sharing one source group
  (max, not sum) so echoing Babel counts once, (b) the guardrail edge-prune,
  which drops cross-family conflation edges before clustering regardless of weight,
  and (c) being **size-aware**: only a small list merges on its own, because it is
  the large ones that are conflated (see ``max_merge_list_size``).
* Correlated sources are de-correlated: KG2 / ROBOKOP / Translator all derive
  equivalence from the SRI Node Normalizer (Babel), so they are **not**
  independent evidence. Evidence within a source group is combined by
  **max**, not sum, so their shared ancestry cannot triple-count and re-import
  Babel's clustering. Independent sources accumulate by sum.

ALL NUMBERS HERE ARE UNTUNED PLACEHOLDERS. They exist so the pipeline runs
end-to-end; the eval harness (``eval/scorer.py``) is what tunes them. Override
via ``config/entity_resolution/weights.yaml``.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import yaml
from pydantic import BaseModel, Field

from kraken.utils.constants import PROJECT_ROOT

# Optional tuning file; if absent, the defaults below apply.
DEFAULT_WEIGHTS_PATH = PROJECT_ROOT / "config" / "entity_resolution" / "weights.yaml"

# Match predicates. exact_match / same_as are full-strength equivalence;
# close_match is weak; broad_match / narrow_match are hierarchical and MUST be
# excluded (including them guarantees parent/child collapse). Predicate values
# are always biolink-prefixed after harmonization, so only prefixed forms appear.
EXACT_MATCH_PREDICATES: frozenset[str] = frozenset({"biolink:exact_match", "biolink:same_as"})
CLOSE_MATCH_PREDICATES: frozenset[str] = frozenset({"biolink:close_match"})
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
            # THE equivalence backbone: the SRI Node Normalizer's cliques (source
            # "nn"), queried live so they're the current, cleanest Babel mapping.
            # Above tau so a clique merges on its own; cross-family conflations are
            # pruned pairwise before clustering, so this only governs SAME-family
            # merges.
            "nn": 0.5,
            # Curated, structurally tight -> reach tau, so they merge on their own.
            "ncbigene": 1.0,
            "refmet": 1.0,
            "lipidmaps": 1.0,
            "umls": 0.8,
            "loinc": 0.6,
            "cdes": 0.6,
            # Aggregators' baked-in equivalent_ids lists. At or above tau, so a list can merge on
            # its own -- but only a SMALL one: these sources are SIZE-AWARE (see
            # max_merge_list_size), which is what actually keeps them safe. Their value is
            # recovering the mappings NN doesn't know, and that gap is large: of kg2's list pairs
            # 64% are unknown to NN (RXCUI/UMLS/RXNORM/CHV), of robokop's 78% (almost all DBSNP).
            # Sharing the "sri_nn_derived" group (max-not-sum) still stops the same Babel
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

    # SIZE-AWARE equivalency. An aggregator's equivalent_ids list is trustworthy when small and
    # conflated when large -- the junk ("kg2 fusing 1300+ Reactome ids into a gene") is a thin tail,
    # not the body. Measured against the live NN, on member pairs NN resolves, agreement by list size:
    #
    #     list size   3-5   6-10  11-20  21-30  31-60  61-100  101+
    #     kg2         94%   96%    95%    73%    37%    19%     8%
    #     robokop    100%   99%   100%    87%    42%    24%    11%
    #
    # Flat through 20, then a cliff -- so rather than distrust every list for the sake of the tail,
    # weight depends on size, in three tiers:
    #
    #   size <= max_merge_list_size[source]   -> equivalency_weights[source]  (at/above tau: merges alone)
    #   size <= max_corroborate_list_size     -> corroborate_weight           (below tau: corroborates)
    #   larger                                -> 0                            (no evidence at all)
    #
    # kg2's lists hold ~95% through 20 and fall to 73% across 21-30; it merges up to 24, reaching a
    # little into that band (the measurement doesn't resolve where within it the drop begins). robokop's
    # degrade more gracefully (87% at 21-30) and its lists never exceed 46, so it merges up to 30.
    # translator-kg-open, unmeasured, gets the conservative 20; it emits no equivalent_ids lists (its
    # equivalent_ids_prop is ""), so for it this only governs aliases.
    #
    # ONLY the sources listed here are size-aware. Everything else keeps its flat weight at any list
    # size -- deliberately, since the curve above was measured on aggregators and says nothing about
    # them. That matters most for "nn": the normalizer's cliques are clean and legitimately large
    # (gene/protein ~30-40), so capping them at 20 would break up the equivalence backbone.
    max_merge_list_size: dict[str, int] = Field(
        default_factory=lambda: {
            "kg2": 24,
            "robokop": 30,
            "translator-kg-open": 20,
        }
    )
    # Mid-size lists (37-87% agreement) are too unreliable to merge on but real enough to
    # corroborate, so they keep the old aggregator weight; beyond 60 (8-24%) they are mostly
    # conflation and contribute nothing. Must stay below tau, or the middle tier would merge alone.
    corroborate_weight: float = 0.15
    max_corroborate_list_size: int = 60
    #
    # An ALIAS -- the original -> canonical pairing a canonicalized aggregator records on an edge
    # (see uncanonicalize.original_alias_pairs) -- is the same kind of assertion as a list entry
    # ("this source says these are one thing"), so it takes the same weight rather than a knob of its
    # own: it is weighted as a list of two. That puts it in the top tier for every source, and means
    # the two can never drift apart.
    ALIAS_LIST_SIZE: ClassVar[int] = 2

    # Per-predicate weights for source match-predicate edges. close_match is a weak
    # "roughly the same" assertion (not true equivalence), so it's kept very low —
    # well below tau, so it only ever corroborates, never merges on its own (and the
    # subclass co-occurrence penalty can drive it lower still).
    #
    # Corroboration here means PARALLEL assertions: several close_match edges on the same pair, each
    # from a different primary knowledge source, are independent claims and should sum -- three at
    # 0.1 reach tau. So match-predicate evidence is grouped per (source, primary knowledge source),
    # not per source (see predicate_group). Grouping by source alone put every kg2 assertion in
    # sri_nn_derived, where they combine by MAX and three parallel edges came to 0.1, not 0.3.
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
            # All Babel/SRI-NN-derived equivalence shares one group so echoing the same
            # Babel assertion (live NN clique + the aggregators' baked-in lists) counts
            # once (max), not summed. "nn" is the live NN cliques; the rest are the
            # aggregators' stored lists.
            "sri_nn_derived": ["nn", "kg2", "robokop", "translator-kg-open"],
        }
    )

    # Equivalency clique handling. Sets up to this size become full cliques; larger
    # sets become a star from the lexically-smallest hub. Now that equivalence comes
    # from the normalizer's CLEAN cliques (not aggregators' junk lists), this is a
    # pure SCALE valve (avoid N^2 edges on a pathological clique), not a distrust
    # mechanism — so keep it above real clique sizes (gene/protein ~30-40, disease
    # ~15-20) so legitimate cliques stay full (robust), and only true outliers star.
    clique_cap: int = 100

    # Name-similarity blocking: skip groups larger than this (likely a generic
    # token), and drop names shorter than this or purely numeric.
    name_group_cap: int = 40
    min_name_length: int = 3

    def model_post_init(self, _context: object) -> None:
        if self.corroborate_weight >= self.tau:
            raise ValueError(
                f"corroborate_weight ({self.corroborate_weight}) must stay below tau ({self.tau}); "
                f"otherwise mid-size equivalency lists would merge on their own"
            )
        # Precompute source -> source-group-id.
        self._source_to_group: dict[str, str] = {}
        for group_id, members in self.source_groups.items():
            for member in members:
                self._source_to_group[member] = group_id

    # ---- lookups ----

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

    def equivalency_weight(self, source: str, list_size: int = ALIAS_LIST_SIZE) -> float:
        """Weight for one clique edge from an equivalency list of ``list_size`` ids.

        Size-aware for the sources in ``max_merge_list_size``; flat for everything else. The default
        size is a pair, i.e. the source's full (top-tier) weight -- which is also what an alias gets.
        """
        base = self.equivalency_weights.get(source, self.default_equivalency_weight)
        max_merge = self.max_merge_list_size.get(source)
        if max_merge is None or list_size <= max_merge:
            return base
        if list_size <= self.max_corroborate_list_size:
            return min(self.corroborate_weight, base)  # never lift a source above its own weight
        return 0.0

    def alias_weight(self, source: str) -> float:
        """Weight for one original -> canonical alias: the same as a two-id equivalency list."""
        return self.equivalency_weight(source, self.ALIAS_LIST_SIZE)

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
