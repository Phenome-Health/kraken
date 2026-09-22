"""Name normalization and name-similarity edges (plan §1).

Name-similarity edges catch entities with no equivalency at all. Rules:

* Exact equality of a **name key** (see ``name_keys``): lowercased, accents and
  possessives dropped, spacing and punctuation removed -- except what changes the
  entity (stereo signs, charges, primes, and the separators between two numbers)
  -- and, outside chemistry, singular forms matched too. Group by the key — no
  blocking, no first-character bucketing (which would silently miss pairs).
* **Primary names only** — synonyms must not enter the graph (NCBI Gene curates
  cleavage products as gene aliases; LOINC attaches ``'Point in time'`` to 97k
  codes).
* Cap group size, keep a stoplist, drop very short and purely numeric names,
  and drop names that are just an identifier (see IDENTIFIER_LIKE_NAME).
* Weighting is handled in ``weights.py`` (well below any single-source
  equivalency), not here.
"""

from __future__ import annotations

import re
import unicodedata
from collections import defaultdict
from collections.abc import Iterable, Iterator

from kraken.entity_resolution.families import ALL_FAMILIES

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


# Names that are really an identifier. A dbSNP rsid names a POSITION, and every ClinGen allele at that position
# is named after it -- 5.1M CAID and 5.0M DBSNP nodes are called "rs10154897" and the like. Matching on those is
# identifier matching, not name matching: it glues every allele at a position into one cluster, which the CAID
# one-id guardrail then has to take apart again (that split is what fills the build log). The allele/position
# relationship is already carried properly, as ``member_of`` edges (see harmonizers/robokop.py).
IDENTIFIER_LIKE_NAME = re.compile(r"^rs\d+$")


def _fold(name: str) -> str:
    """Lowercase and strip combining accents."""
    text = unicodedata.normalize("NFKD", name)
    return "".join(ch for ch in text if not unicodedata.combining(ch)).lower()


def normalize_name(name: str | None) -> str:
    """Lowercase, strip accents, collapse punctuation/whitespace to single spaces.

    The WORDS of a name: what ``is_droppable`` judges, and the matching key for the branches whose names differ
    only in punctuation (see ``name_keys``)."""
    if not name:
        return ""
    return _PUNCT_WS.sub(" ", _fold(name)).strip()


# Families whose names keep their spacing and punctuation in the matching key. Distinct NCBI taxa differ only by it
# ("Burkholderia sp. S2" / "Burkholderia sp. S-2"), and gene, protein and variant names are symbols and HGVS strings,
# where it is part of the identifier.
SPACED_NAME_FAMILIES: frozenset[str] = frozenset({"organism", "gene_protein", "genomic_variant"})
# Families whose names are never singularized: chemistry uses the plural for a CLASS ("uridines" is not "uridine",
# "Resorcinols" is not resorcinol). Names of a wildcard node aren't either: it may be a chemical.
UNSINGULARIZED_FAMILIES: frozenset[str] = frozenset({"chemical"})

# Kept in the key because they change the entity: a stereo sign or charge in parentheses ("(+)-camphor", "(2-)"),
# and a prime ("3'-glucoside" is not "3-glucoside"). "(+/-)", "(+-)" and "(±)" all mean racemic.
_SIGN = re.compile(r"\((\d*[+-]|\+/-|\+-|±)\)")
_POSSESSIVE = re.compile(r"(?<=[^\W\d_])'s\b")
_PRIMES = str.maketrans({"’": "'", "′": "'", "`": "'", "″": "''"})
_TOKEN = re.compile(r"\((?:\d*[+-]|±)\)|'|[^\W_]+", re.UNICODE)
# Punctuation between two numbers that says how they relate: a decimal point ("0.1" is not the range "0-1"), and a
# lipid's chains with their sn-positions known ("16:0/18:1") or not ("16:0_18:1"). Anything else there is a plain
# break between two locants, written many ways ("14,15-", "14(15)-", "1-1-1-", a reaction's "O2 => (3S)" or
# "O2 <=> (3S)"), and keyed as ",". What matters is that it is there at all: "1,2-" is not "12-".
_NUMBER_SEPARATORS = frozenset("./:_")
# Anything in a raw name that the compact key keeps and the spaced words lose.
_DISTINGUISHING = re.compile(r"\((?:\d*[+-]|\+/-|\+-|±)\)|['’′`″]|\d[./:_]\d")
# Word endings that are not plurals: "status", "fibrosis", "venous", "glass".
_NOT_PLURAL = ("ss", "us", "is", "ous")


def _singular(word: str) -> str:
    if len(word) <= 3 or any(ch.isdigit() for ch in word) or word.endswith(_NOT_PLURAL):
        return word
    if word.endswith("ies") and len(word) > 4:
        return word[:-3] + "y"
    return word[:-1] if word.endswith("s") else word


def _compact(name: str, *, singularize: bool) -> str:
    text = _POSSESSIVE.sub("", _fold(name).translate(_PRIMES))
    text = _SIGN.sub(lambda m: "(±)" if m.group(1) in ("+/-", "+-") else m.group(0), text)
    out: list[str] = []
    signs: list[str] = []
    previous, previous_end = "", 0
    for match in _TOKEN.finditer(text):
        token = match.group()
        if token.startswith("("):
            signs.append(token)  # position-free: "(+)-limonene" is "LIMONENE, (+)-"
            continue
        if previous[-1:].isdigit() and token[:1].isdigit():
            between = [ch for ch in text[previous_end : match.start()] if ch in _NUMBER_SEPARATORS]
            out.append(between[0] if between else ",")
        if singularize and token[:1].isalpha():
            token = _singular(token)
        out.append(token)
        previous, previous_end = token, match.end()
    return "".join(out) + "".join(sorted(signs))


def name_keys(name: str | None, branches: frozenset[str] = ALL_FAMILIES) -> frozenset[str]:
    """The keys a name matches on, given its node's branch families; two names match if they share one.

    The main key drops spacing, hyphens and other punctuation ("LY-2940094" = "LY2940094", "AM 404" = "AM404"), and
    a possessive ("Parkinson's disease" = "Parkinson disease"). What distinguishes one entity from another stays: a
    stereo sign or charge in parentheses (wherever it is written), a prime, and the separator between two numbers (so
    "1,2-" is not "12-", and the lipid "PC 14:0/20:0", sn-positions known, is not "PC 14:0_20:0").

    Two more keys, so no pair of identical names can end up keyed apart:

    * a typed node outside chemistry gets a SINGULAR key too ("Jejunal Neoplasms" -> also "jejunal neoplasm"), which
      meets the singular name on any node. A chemical plural is a class, and gets none; nor does an untyped node's,
      which may be a chemical.
    * organisms, genes, proteins and variants key on the spaced words of ``normalize_name`` only, and a node with no
      type (a wildcard) keys on those as well, so an untyped "Mobiluncus sp" still meets NCBITaxon's -- unless its
      name has a sign, charge, prime, decimal point or lipid separator, which the spaced words would lose.
    """
    if not name:
        return frozenset()
    spaced = normalize_name(name)
    if branches is not ALL_FAMILIES and branches & SPACED_NAME_FAMILIES:
        return frozenset({spaced})
    keys = {_compact(name, singularize=False)}
    if branches is ALL_FAMILIES:
        if not _DISTINGUISHING.search(name):
            keys.add(spaced)
    elif not branches & UNSINGULARIZED_FAMILIES:
        keys.add(_compact(name, singularize=True))
    return frozenset(keys)


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
    if IDENTIFIER_LIKE_NAME.match(normalized):
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
