"""
Labels for callers, added to a built graph's nodes. They change no id, merge, edge or other field.

pathway_template_family: PathWhiz makes one pathway per lipid species from a template, like "Cardiolipin Biosynthesis
CL(10:0/10:0/10:0/28:0) pathway". 93,208 such copies in 9 families share most of their edges, so they crowd out
everything else in a result. The label names the family, the copy's own lipid part and the family's size, so a caller
can group copies before cutting results and still see which lipid each one is about. It is read from the name;
PathWhiz itself doesn't mark them.

same_name_chemical_candidates: a node with only vocabulary ids (a lab test, an answer list, a concept typed wrong)
and only mapping, hierarchy or mention edges, named exactly like a chemical that has biological edges: "Ozone"
LOINC:105571-4 next to "Ozone" HMDB:HMDB0035409. The label lists those chemicals as candidates a caller can prefer. It
never says the two are the same thing.

Run on a built graph:  python -m kraken.post_processing.node_labels NODES_JSONL EDGES_JSONL
"""

import argparse
import json
import logging
import os
import re
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from kraken.utils.constants import (
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_SUBJECT,
    KRAKEN_SOURCE_ID,
    NODE_ATTRIBUTES,
    NODE_CATEGORIES,
    NODE_EQUIVALENT_IDS,
    NODE_ID,
    NODE_NAME,
)

FAMILY_ATTRIBUTE = "pathway_template_family"
SAME_NAME_ATTRIBUTE = "same_name_chemical_candidates"
LABELS = (FAMILY_ATTRIBUTE, SAME_NAME_ATTRIBUTE)

TEMPLATE_PATHWAY_PREFIXES = ("PathWhiz:", "SMPDB:")
# Where a lipid's shorthand starts in a pathway name; the chains themselves aren't parsed ("CL(i-19:0/10:0(3-OH)/...").
LIPID_START = re.compile(r"\b(?:TG|DG|MG|PC|PE|PS|PA|PG|PI|CL|BMP|SM|Cer|CDP-DG|LPC|LPE|LysoPC|LysoPE|PIP|PIP2)\s*\(")

VOCABULARY_PREFIXES = frozenset(
    {"LOINC", "UMLS", "CHV", "MEDDRA", "SNOMEDCT", "NCIT", "MESH", "RXCUI", "NDDF", "VANDF", "ATC"}
)
CHEMICAL_PREFIXES = frozenset(
    {
        "INCHIKEY", "CHEBI", "PUBCHEM.COMPOUND", "HMDB", "KEGG", "KEGG.COMPOUND", "RM", "LM", "CAS", "CHEMBL.COMPOUND",
        "UNII", "DRUGBANK", "SMILES",
    }
)  # fmt: skip


def _biolink(*names: str) -> frozenset[str]:
    return frozenset(f"biolink:{name}" for name in names)


PATHWAY_CATEGORY = "biolink:Pathway"
CHEMICAL_CATEGORIES = _biolink(
    "SmallMolecule", "ChemicalEntity", "Drug", "MolecularMixture", "ChemicalMixture", "NucleicAcidEntity"
)
# A chemical "has biology" when it has one of these edges; a vocabulary node keeps its label only if all its edges are
# the other kind.
BIOLOGICAL_PREDICATES = _biolink(
    "affects", "interacts_with", "physically_interacts_with", "regulates", "has_participant", "actively_involved_in",
    "participates_in", "produces", "consumes", "has_input", "has_output", "catalyzes", "enables", "occurs_in",
)  # fmt: skip
VOCABULARY_PREDICATES = _biolink(
    "subclass_of", "superclass_of", "same_as", "exact_match", "close_match", "broad_match", "narrow_match",
    "related_to", "mentions",
)  # fmt: skip


def _prefix(curie: str) -> str:
    return curie.split(":", 1)[0]


def _ids(node: dict[str, Any]) -> set[str]:
    return {node[NODE_ID], *node.get(NODE_EQUIVALENT_IDS, ())}


def template_family(node: dict[str, Any]) -> tuple[str, str] | None:
    """(family, this copy's lipid part) of a PathWhiz or SMPDB pathway made from a lipid template, else None."""
    if PATHWAY_CATEGORY not in node.get(NODE_CATEGORIES, ()):
        return None
    if not any(curie.startswith(TEMPLATE_PATHWAY_PREFIXES) for curie in _ids(node)):
        return None
    name = node.get(NODE_NAME) or ""
    match = LIPID_START.search(name)
    family = name[: match.start()].strip() if match else ""
    return (family, name[match.start() :].removesuffix(" pathway")) if family else None


def is_vocabulary_only(node: dict[str, Any]) -> bool:
    return bool(node.get(NODE_NAME)) and all(_prefix(curie) in VOCABULARY_PREFIXES for curie in _ids(node))


def is_chemical(node: dict[str, Any]) -> bool:
    return not CHEMICAL_CATEGORIES.isdisjoint(node.get(NODE_CATEGORIES, ())) and any(
        _prefix(curie) in CHEMICAL_PREFIXES for curie in _ids(node)
    )


def _nodes(path: Path) -> Iterator[dict[str, Any]]:
    with open(path) as handle:
        for line in handle:
            yield json.loads(line)


def _edge_ends(path: Path) -> Iterator[tuple[str, str, str]]:
    with open(path) as handle:
        for line in handle:
            edge = json.loads(line)
            yield edge[EDGE_SUBJECT], edge[EDGE_PREDICATE], edge[EDGE_OBJECT]


def find_labels(nodes_path: Path, edges_path: Path) -> tuple[Counter, dict[str, list[str]]]:
    """Family sizes, and each labelled vocabulary node's chemical candidates. Two node passes and one edge pass; memory
    holds the vocabulary-only names (about 2.2M) for one pass, then only the few thousand that match a chemical."""
    family_size: Counter = Counter()
    vocabulary_by_name: dict[str, list[str]] = defaultdict(list)
    for node in _nodes(nodes_path):
        family = template_family(node)
        if family:
            family_size[family[0]] += 1
        elif is_vocabulary_only(node):
            vocabulary_by_name[node[NODE_NAME]].append(node[NODE_ID])

    chemicals_by_name: dict[str, list[str]] = defaultdict(list)
    for node in _nodes(nodes_path):
        if node.get(NODE_NAME) in vocabulary_by_name and is_chemical(node):
            chemicals_by_name[node[NODE_NAME]].append(node[NODE_ID])
    vocabulary = {curie: name for name in chemicals_by_name for curie in vocabulary_by_name[name]}
    chemicals = {curie for ids in chemicals_by_name.values() for curie in ids}
    del vocabulary_by_name

    biological: set[str] = set()
    more_than_vocabulary: set[str] = set()
    for subject, predicate, object_ in _edge_ends(edges_path):
        for curie in (subject, object_):
            if curie in chemicals and predicate in BIOLOGICAL_PREDICATES:
                biological.add(curie)
            if curie in vocabulary and predicate not in VOCABULARY_PREDICATES:
                more_than_vocabulary.add(curie)

    candidates = {}
    for curie, name in vocabulary.items():
        chemical_ids = sorted(c for c in chemicals_by_name[name] if c in biological and c != curie)
        if chemical_ids and curie not in more_than_vocabulary:
            candidates[curie] = chemical_ids
    return family_size, candidates


def labelled(node: dict[str, Any], family_size: Counter, candidates: dict[str, list[str]]) -> dict[str, Any] | None:
    """The node with its labels set (and stale ones removed), or None when nothing about it changes."""
    labels: dict[str, Any] = {}
    family = template_family(node)
    if family:
        labels[FAMILY_ATTRIBUTE] = {"name": family[0], "variant": family[1], "members": family_size[family[0]]}
    if node[NODE_ID] in candidates:
        labels[SAME_NAME_ATTRIBUTE] = candidates[node[NODE_ID]]
    attributes = node.get(NODE_ATTRIBUTES) or {}
    ours = attributes.get(KRAKEN_SOURCE_ID) or {}
    kept = {key: value for key, value in ours.items() if key not in LABELS}
    if {**kept, **labels} == ours:
        return None
    attributes = {**attributes, KRAKEN_SOURCE_ID: {**kept, **labels}}
    if not attributes[KRAKEN_SOURCE_ID]:
        del attributes[KRAKEN_SOURCE_ID]
    node = {key: value for key, value in node.items() if key != NODE_ATTRIBUTES}
    return {**node, NODE_ATTRIBUTES: attributes} if attributes else node


def label_nodes(nodes_path: Path, edges_path: Path) -> Counter:
    """Rewrite the nodes file with the labels, through a new temporary file beside it. Edges are only read."""
    family_size, candidates = find_labels(nodes_path, edges_path)
    counts: Counter = Counter()
    descriptor, name = tempfile.mkstemp(dir=nodes_path.parent, prefix=f"{nodes_path.name}.", suffix=".labelling")
    temporary = Path(name)
    try:
        with open(nodes_path) as nodes_in, os.fdopen(descriptor, "w") as nodes_out:
            for line in nodes_in:
                counts["nodes"] += 1
                node = labelled(json.loads(line), family_size, candidates)
                if node is None:
                    nodes_out.write(line)
                    continue
                counts["changed"] += 1
                nodes_out.write(json.dumps(node, ensure_ascii=False) + "\n")
            nodes_out.flush()
            os.fsync(nodes_out.fileno())  # the new file must be on disk before it replaces the build's
        os.chmod(temporary, nodes_path.stat().st_mode)  # mkstemp makes it private
        os.replace(temporary, nodes_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    counts[FAMILY_ATTRIBUTE] = sum(family_size.values())
    counts[SAME_NAME_ATTRIBUTE] = len(candidates)
    logging.info(
        "node labels: %d pathway template copies in %d families, %d vocabulary nodes with same-name chemical "
        "candidates; %d of %d nodes rewritten",
        counts[FAMILY_ATTRIBUTE],
        len(family_size),
        counts[SAME_NAME_ATTRIBUTE],
        counts["changed"],
        counts["nodes"],
    )
    return counts


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("nodes", type=Path)
    parser.add_argument("edges", type=Path)
    arguments = parser.parse_args()
    label_nodes(arguments.nodes, arguments.edges)
