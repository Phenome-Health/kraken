"""Per-source / per-predicate evidence weights for the match graph.

Design (see ``docs/entity_resolution_plan.md`` §1):

* Weights are **per source** and **accumulate across agreeing sources**.
* Discriminative power comes from sources *disagreeing*, so weight design
  matters more than gamma. The key calibration is relative to gamma (the
  clustering resolution): a source whose weight exceeds gamma can merge a pair
  on its own; a source below gamma can only *corroborate*. Over-conflating
  sources (the aggregators) therefore sit **below** gamma.
* Correlated sources are de-correlated: KG2 / ROBOKOP / Translator all derive
  equivalence from the SRI Node Normalizer (Babel), so they are **not**
  independent evidence. Evidence within a correlation group is combined by
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

    # Clustering resolution. A and B joined by total weight w stay separate iff
    # gamma * |A| * |B| > w. Kept here so weights can be reasoned about relative
    # to it (the file that documents the calibration owns the number).
    gamma: float = 0.5

    # Minimum accumulated pair weight to keep an edge at all (pre-filter). Runs
    # before connected components, so it shapes the component structure.
    tau: float = 0.3

    # Per-source weight for an equivalency-list assertion (one clique edge).
    # Sources absent from this map use ``default_equivalency_weight``.
    equivalency_weights: dict[str, float] = Field(
        default_factory=lambda: {
            # Curated, structurally tight -> can merge on their own (> gamma).
            "ncbigene": 1.0,
            "refmet": 1.0,
            "lipidmaps": 1.0,
            "umls": 0.8,
            "loinc": 0.6,
            "cdes": 0.6,
            # Aggregators / SRI-NN-derived -> corroboration only (< gamma).
            "kg2": 0.3,
            "robokop": 0.3,
            "translator-kg-open": 0.3,
            "microbiome-kg": 0.3,
            "multiomics-kg": 0.3,
        }
    )
    default_equivalency_weight: float = 0.4

    # Per-predicate weights for source match-predicate edges.
    exact_match_weight: float = 1.0
    close_match_weight: float = 0.2

    # Name-similarity edges: ABOVE gamma, so two nodes linked ONLY by a shared
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

    # Correlation groups: sources listed together contribute by MAX, not sum.
    # Each source maps to a group id; sources not listed are their own group.
    correlation_groups: dict[str, list[str]] = Field(
        default_factory=lambda: {
            "sri_nn_derived": ["kg2", "robokop", "translator-kg-open", "microbiome-kg", "multiomics-kg"],
        }
    )

    # Equivalency clique handling. Sets up to this size become full cliques;
    # larger sets become stars (fragile on purpose — big "equivalent" lists are
    # usually junk). See plan §1.
    clique_cap: int = 20

    # Name-similarity blocking: skip groups larger than this (likely a generic
    # token), and drop names shorter than this or purely numeric.
    name_group_cap: int = 40
    min_name_length: int = 3

    def model_post_init(self, _context: object) -> None:
        # Precompute source -> correlation-group-id.
        self._source_to_group: dict[str, str] = {}
        for group_id, members in self.correlation_groups.items():
            for member in members:
                self._source_to_group[member] = group_id

    # ---- lookups ----

    def correlation_group(self, source: str) -> str:
        """Group id used for de-correlation; independent sources get their own."""
        return self._source_to_group.get(source, f"src:{source}")

    def equivalency_weight(self, source: str) -> float:
        return self.equivalency_weights.get(source, self.default_equivalency_weight)

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
