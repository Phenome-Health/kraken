"""
Name overrides for vocabularies whose names don't say what kind of thing the node is.

Some vocabularies name a node after the thing it's *about* rather than what it *is*: KEGG's Parkinson disease
pathway is called "Parkinson disease", LOINC's answer code LA27533-1 is called "Parkinson's disease", and a
Rhea reaction is named for its enzyme. Alongside the real disease node these are indistinguishable in a result
list, so a word is appended to say what they are -- "Parkinson disease pathway", "Parkinson's disease [answer]".

Rules are keyed on (id prefix, category) rather than category alone, because the right treatment differs by
vocabulary. GO's MolecularActivity names are already 82% self-identifying ("protein kinase activity"), so a
blanket per-category rule would produce "...activity activity".

Two properties every rule must keep:

  - APPEND-ONLY. The original name stays a substring of the new one, so text search still finds it. That is
    what makes it safe NOT to keep the bare name as a synonym: an exact-match synonym would rank the answer
    code alongside the real disease, which is the confusion this exists to remove. A rule that reordered or
    rewrote the name would break that and would need the synonym back.
  - IDEMPOTENT. Skipped when the name already contains the keyword, so re-running is safe and vocabularies
    that self-identify some of the time aren't double-suffixed.

The pre-override name is kept in attributes under ORIGINAL_NAME_ATTRIBUTE.
"""

from dataclasses import dataclass

# Attribute holding the name as the source gave it, recorded only on nodes that were actually renamed.
ORIGINAL_NAME_ATTRIBUTE = "original_name"


@dataclass(frozen=True)
class NameOverride:
    """suffix -- appended after a space ("Glycolysis" -> "Glycolysis pathway").
    already_says_it -- substrings (matched case-insensitively) that mean the name already conveys what the
    node is, so it's left alone. Usually just the suffix word, but it also carries shapes that are
    self-evident without naming the type: a Rhea name containing " + " or "=" is visibly a chemical equation,
    and 98% of them are, so only the handful named after their enzyme get suffixed."""

    suffix: str
    already_says_it: tuple[str, ...]


# (id prefix, category) -> override. The id prefix is matched with startswith, so a vocabulary that splits its
# namespace by local-id pattern can be targeted precisely (LOINC's LA answers vs its LP parts and LL lists).
# Node counts are from KRAKEN 2.1.0; the "self-identifying" share of each is near zero unless noted.
NAME_OVERRIDES: dict[tuple[str, str], NameOverride] = {
    # --- Pathways named after the disease/process they describe ---
    ("SMPDB:", "biolink:Pathway"): NameOverride("pathway", ("pathway",)),  # 62,473
    ("PathWhiz:", "biolink:Pathway"): NameOverride("pathway", ("pathway",)),  # 32,293
    ("REACT:", "biolink:Pathway"): NameOverride("pathway", ("pathway",)),  # 20,955
    ("KEGG:", "biolink:Pathway"): NameOverride("pathway", ("pathway",)),  # the issue's example, KEGG:05012
    ("GO:", "biolink:Pathway"): NameOverride("pathway", ("pathway",)),  # 1,163
    # --- Gene families named identically to one of their members ---
    # The largest source of confusion by volume: 18,918 PANTHER families share a name with a real gene or
    # protein node ("carbonic anhydrase", "caspase", "monoamine oxidase").
    ("PANTHER.FAMILY:", "biolink:GeneFamily"): NameOverride("family", ("family",)),  # 26,139
    ("HGNC.FAMILY:", "biolink:GeneFamily"): NameOverride("family", ("family",)),  # 1,706
    # --- Reactions named after the enzyme that catalyses them ---
    # Only 9 of Rhea's 520 nodes are named this way ("creatininase", "xylulokinase"); the other 98% are named
    # by their chemical equation, which already reads as a reaction, so " + " and "=" suppress the suffix.
    ("RHEA:", "biolink:MolecularActivity"): NameOverride("reaction", ("reaction", " + ", "=")),
    # --- Answer codes carrying the verbatim answer text ---
    # Bracketed rather than suffixed: these are information artifacts, and "Parkinson's disease answer" would
    # read as a kind of answer rather than as a code whose text happens to be a disease name.
    ("LOINC:LA", "biolink:InformationContentEntity"): NameOverride("[answer]", ("[answer]",)),  # 22,146
}

# Indexed by category so the common case -- a node whose category has no rule at all -- costs one dict miss
# rather than a scan of every id prefix in the table.
_OVERRIDES_BY_CATEGORY: dict[str, tuple[tuple[str, NameOverride], ...]] = {}
for (_id_prefix, _category), _override in NAME_OVERRIDES.items():
    _OVERRIDES_BY_CATEGORY.setdefault(_category, ())
    _OVERRIDES_BY_CATEGORY[_category] += ((_id_prefix, _override),)


def find_name_override(curie: str, categories: list[str]) -> NameOverride | None:
    """The override for this node, or None if no rule applies."""
    for category in categories:
        for id_prefix, override in _OVERRIDES_BY_CATEGORY.get(category, ()):
            if curie.startswith(id_prefix):
                return override
    return None


def apply_name_override(curie: str, categories: list[str], name: str) -> str:
    """The node's name with its type appended, or the name unchanged when no rule applies or the name already
    says what it is."""
    override = find_name_override(curie, categories)
    if override is None:
        return name
    lowered = name.lower()
    if any(marker.lower() in lowered for marker in override.already_says_it):
        return name
    return f"{name} {override.suffix}"
