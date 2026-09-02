"""Name normalization and name-similarity edges (plan §1).

Name-similarity edges catch entities with no equivalency at all. Rules:

* Exact equality after normalization (lowercase, strip punctuation/whitespace).
  Group by the normalized string — no blocking, no first-character bucketing
  (which would silently miss pairs).
* **Primary names only** — synonyms must not enter the graph (NCBI Gene curates
  cleavage products as gene aliases; LOINC attaches ``'Point in time'`` to 97k
  codes).
* Cap group size, keep a stoplist, drop very short and purely numeric names.
* Weighting is handled in ``weights.py`` (well below any single-source
  equivalency), not here.
"""

from __future__ import annotations

import re
import unicodedata
from collections import defaultdict
from collections.abc import Iterable, Iterator

_PUNCT_WS = re.compile(r"[^\w]+", re.UNICODE)

# Generic tokens that collide across unrelated entities. Extend as the eval
# surfaces offenders (plan calls for a stoplist).
DEFAULT_STOPLIST: frozenset[str] = frozenset(
    {
        "point in time",
        "unknown",
        "other",
        "none",
        "not applicable",
        "normal",
        "abnormal",
        "present",
        "absent",
        "positive",
        "negative",
    }
)


def normalize_name(name: str | None) -> str:
    """Lowercase, strip accents, collapse punctuation/whitespace to single spaces."""
    if not name:
        return ""
    text = unicodedata.normalize("NFKD", name)
    text = "".join(ch for ch in text if not unicodedata.combining(ch))
    text = _PUNCT_WS.sub(" ", text.lower()).strip()
    return text


def is_droppable(normalized: str, *, min_length: int = 3, stoplist: frozenset[str] = DEFAULT_STOPLIST) -> bool:
    """True if a normalized name is too weak to key a similarity group."""
    if not normalized:
        return True
    if len(normalized) < min_length:
        return True
    if normalized in stoplist:
        return True
    # purely numeric (digits and spaces only)
    if all(ch.isdigit() or ch.isspace() for ch in normalized):
        return True
    return False


def group_by_normalized_name(
    id_name_pairs: Iterable[tuple[str, str | None]],
    *,
    min_length: int = 3,
    stoplist: frozenset[str] = DEFAULT_STOPLIST,
) -> dict[str, list[str]]:
    """Group CURIEs by normalized primary name, dropping weak names.

    ``id_name_pairs`` yields ``(curie, primary_name)`` — pass primary names only.
    """
    groups: dict[str, set[str]] = defaultdict(set)
    for curie, name in id_name_pairs:
        norm = normalize_name(name)
        if is_droppable(norm, min_length=min_length, stoplist=stoplist):
            continue
        groups[norm].add(curie)
    return {norm: sorted(ids) for norm, ids in groups.items() if len(ids) > 1}


def name_similarity_edges(
    groups: dict[str, list[str]],
    *,
    group_cap: int = 40,
) -> Iterator[tuple[str, str]]:
    """Yield unordered CURIE pairs (a < b) that share a normalized name.

    Groups larger than ``group_cap`` are skipped (a generic token, likely junk).
    """
    for _norm, ids in groups.items():
        if len(ids) > group_cap:
            continue
        n = len(ids)
        for i in range(n):
            for j in range(i + 1, n):
                a, b = ids[i], ids[j]
                yield (a, b) if a < b else (b, a)
