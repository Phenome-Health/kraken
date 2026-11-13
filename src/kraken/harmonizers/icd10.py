# icd10.py
"""ICD-10 -- WHO's own classification, the diagnosis codes of death records and of most patient records outside the
US -- beside ICD-10-CM, its US modification (see icd10cm).

Input is a directory (config ``input_file``) holding the code file of one WHO edition. Two things come from the
harmonized output of the build's other sources: Babel's edges, and the nodes of the icd10cm source.

  * every code Babel does not already place is a node with a ``biolink:subclass_of`` edge to its parent. The tree is
    in the code itself, as for ICD-10-CM.
  * a code Babel has in a clique is left as it is: no node, no name, no edge from it or to it. Babel holds 1,841 of
    WHO's codes in disease cliques, as cross-references that are often not the disease (``Q93.5`` "Other deletions of
    part of a chromosome" on distal 10q deletion syndrome), and that is where they have resolved since before this
    source. Named here, such a code would bring the ids of that name in after it, its ICD-10-CM twin first. An edge
    to its parent would hang WHO's tree on the disease (syphilis under "Early syphilis", through ``A51.0``), and an
    edge from a child would hang the child on it (``M30.8`` under polyarteritis nodosa, where Babel has ``M30``
    "Polyarteritis nodosa and related conditions"). So the 1,161 children of such codes have no parent edge.
  * a code is the ICD-10-CM code with the same code only where both have one title, British spelling aside: a
    ``biolink:same_as`` edge. Entity resolution does not weigh it as evidence: it puts the code into the node of that
    ICD-10-CM code after clustering, as one more id of that code's type (see ``entity_resolution/icd10_twins``). The
    code alone is not enough: ``E11.0`` is "Type 2 diabetes mellitus: With coma" here and "... with hyperosmolarity"
    there. Of the 10,535 codes both classifications have (2019 edition, FY2027), 7,019 have one title, 609 of them
    once British spelling is folded. Such a code carries its name as one entity resolution does not match on:
    matched on, the name is evidence that draws the ids of that name, its twin among them, into one node with the
    code, and moving the code afterwards leaves them there. For 300 of these twins that is an id their ICD-10-CM
    code is kept apart from, because it merges through Mondo (see ``icd10cm._code_node``). The twin's node is named
    anyway, and a code that does not join it keeps this name.
  * a code titled by the drug or the procedure alone carries that title as a name not matched on either (see
    ``AGENT_TITLED_FROM``).

Nothing links a code to a disease: WHO's file says nothing of that. A code with a name to match on joins what is
named exactly like it, as any id does.

The ICD-10-CM twin at the other end of a ``same_as`` edge gets a stub node, as every harmonized edge needs both its
nodes: an id with no name and no type of its own, so it says nothing about what the icd10cm source describes.

What was read of the other sources is recorded beside the output (``READS_FILENAME``), and the build refuses to
integrate output made from other facts than its Babel and icd10cm output hold now (``output_is_current``): a run
often harmonizes one source only.
"""

import hashlib
import json
import logging
import re
import unicodedata
from collections import Counter
from pathlib import Path
from typing import NamedTuple

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.icd10cm import (
    CATEGORY_BY_CHAPTER,
    DISEASE_CATEGORY,
    SHORTEST_CODE_LENGTH,
    SUBCLASS_OF_PREDICATE,
    cm_curie,
    dotted,
    harmonized_titles,
    parent_code,
    single_file,
)
from kraken.utils.constants import (
    AUTOMATED_AGENT,
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_SUBJECT,
    KNOWLEDGE_ASSERTION,
    KRAKEN_SOURCE_ID,
    MANUAL_AGENT,
    SAME_AS_PREDICATE,
)
from kraken.utils.kg_io import save_to_jsonl

CODES_FILE_PATTERN = "icd10*syst_codes.txt"
ICD10_PREFIX = "ICD10"
ICD10CM_INFORES = "infores:icd10cm"
LEAF_ATTRIBUTE = "leaf"
STUB_CATEGORY = "biolink:NamedThing"
READS_FILENAME = "reads.json"

# WHO's code file: 17 fields split by ";". (WHO's metainfo page describes an older layout.)
FIELD_SEPARATOR = ";"
FIELD_COUNT = 17
LEVEL_FIELD, LEAF_FIELD, BLOCK_FIELD, CODE_FIELD, TITLE_FIELD, CATEGORY_TITLE_FIELD = 0, 1, 4, 7, 8, 9
LEAF_FLAGS = {"T": True, "N": False}
CODE = re.compile(r"[A-Z]\d{2,4}")

# British spellings WHO uses and ICD-10-CM does not, each with the titles of one code it alone makes equal (2019
# edition against FY2027): anaemia / anemia 379, oedema / edema 90, labour / labor 94, goitre / goiter 17,
# counselling / counseling 11. A fold can afford to be greedy: it only decides between the two titles of ONE code.
BRITISH_SPELLING: tuple[tuple[re.Pattern, str], ...] = (
    (re.compile(r"ae"), "e"),
    (re.compile(r"oe"), "e"),
    (re.compile(r"(?<=[a-z]{3})our"), "or"),
    (re.compile(r"(?<=[a-z]{3})re\b"), "er"),
    (re.compile(r"ll(?=ed\b|ing\b)"), "l"),
)


# Two blocks title a code by the drug or the procedure alone, though the code is the harm that came of it. Block,
# named by its first category -> the code length from which that holds:
#   * Y40-Y59, "Drugs, medicaments and biological substances causing adverse effects in therapeutic use": every
#     code. ``Y44.2`` is "Anticoagulants".
#   * Y83-Y84, procedures "as the cause of abnormal reaction of the patient": the subdivisions. ``Y84.0`` is
#     "Cardiac catheterization"; the titles of ``Y83`` and ``Y84`` themselves say it in full.
# Matched on, such a title is the drug class or the procedure: in the 2.4.0-dev5 build ``Y84.0`` joined MeSH's and
# NCIT's cardiac catheterization, and ``Y44.2`` put LOINC's "Anticoagulant" in an exact tie between it and CHEBI's
# anticoagulant, which cost CHEBI's node that id.
AGENT_TITLED_FROM: dict[str, int] = {"Y40": 3, "Y83": 4}


class CodeEntry(NamedTuple):
    title: str
    category_title: str  # the title of the code's 3-character category
    leaf: bool
    block: str  # the first category of the code's block


def read_codes_file(path: Path) -> dict[str, CodeEntry]:
    """``code without its dot -> entry``. A code's level is its length; its eighth field is the code with no dot,
    dagger or asterisk."""
    entries: dict[str, CodeEntry] = {}
    with open(path, encoding="utf-8") as f:
        for number, line in enumerate(f, start=1):
            fields = line.rstrip("\r\n").split(FIELD_SEPARATOR)
            if len(fields) != FIELD_COUNT:
                raise ValueError(f"{path}, line {number}: not a code file line: {line!r}")
            level, leaf, code, title = (fields[i] for i in (LEVEL_FIELD, LEAF_FIELD, CODE_FIELD, TITLE_FIELD))
            if not CODE.fullmatch(code) or str(len(code)) != level or leaf not in LEAF_FLAGS or not title:
                raise ValueError(f"{path}, line {number}: not a code file line: {line!r}")
            if code in entries:
                raise ValueError(f"{path}, line {number}: code {code} is listed twice")
            entries[code] = CodeEntry(title, fields[CATEGORY_TITLE_FIELD], LEAF_FLAGS[leaf], fields[BLOCK_FIELD])
    if not entries:
        raise ValueError(f"{path} lists no code")
    return entries


def code_names(entries: dict[str, CodeEntry]) -> dict[str, str]:
    """``code -> its name``: its title, under its category's where the title alone names several codes. WHO prints a
    subdivision under its category, so 23 titles stand for 75 codes: "By spouse or partner" is ``Y06.0`` under
    "Neglect and abandonment" and ``Y07.0`` under "Other maltreatment". Left so, entity resolution would merge them
    on the name."""
    codes_with_title = Counter(entry.title for entry in entries.values())
    names = {
        code: (
            entry.title
            if codes_with_title[entry.title] == 1 or len(code) == SHORTEST_CODE_LENGTH
            else f"{entry.category_title}: {entry.title}"
        )
        for code, entry in entries.items()
    }
    if shared := sorted(name for name, count in Counter(names.values()).items() if count > 1):
        raise ValueError(f"WHO codes still share a name under their category's title: {shared}")
    return names


def title_key(title: str) -> str:
    """A title as it compares between WHO's list and ICD-10-CM's: no case, accent, punctuation or British spelling
    ("Guillain-Barré syndrome" is "Guillain-Barre syndrome" there, "Anaemia, unspecified" is "Anemia, unspecified")."""
    unaccented = "".join(ch for ch in unicodedata.normalize("NFKD", title) if not unicodedata.combining(ch))
    key = " ".join(re.sub(r"[^a-z0-9]+", " ", unaccented.lower()).split())
    for british, american in BRITISH_SPELLING:
        key = british.sub(american, key)
    return key


def titled_by_agent(code: str, entry: CodeEntry) -> bool:
    shortest = AGENT_TITLED_FROM.get(entry.block)
    return shortest is not None and len(code) >= shortest


def who_curie(code: str) -> str:
    return f"{ICD10_PREFIX}:{dotted(code)}"


def codes_babel_places(babel_edges: Path | None) -> set[str]:
    """The ``ICD10:`` ids Babel has in a clique: the ends of its ``same_as`` edges (see harmonizers/babel.py). No id
    where the build has no Babel (``babel_edges`` is None). A Babel that is in the build but not harmonized yet is an
    error: every code would come out as if nothing placed it."""
    placed: set[str] = set()
    if babel_edges is None:
        return placed
    if not babel_edges.exists():
        raise FileNotFoundError(f"{babel_edges}: Babel must be harmonized before WHO's ICD-10")
    marker = f'"{ICD10_PREFIX}:'
    with open(babel_edges, encoding="utf-8") as f:
        for line in f:
            if marker not in line:
                continue
            edge = json.loads(line)
            if edge[EDGE_PREDICATE] == SAME_AS_PREDICATE:
                ends = (edge[EDGE_SUBJECT], edge[EDGE_OBJECT])
                placed.update(end for end in ends if end.startswith(f"{ICD10_PREFIX}:"))
    return placed


def icd10cm_titles(icd10cm_nodes: Path | None) -> dict[str, str]:
    """``code -> title`` of the build's ICD-10-CM codes, from the icd10cm source's harmonized nodes: the codes this
    build will hold. No code where the build has no icd10cm (``icd10cm_nodes`` is None); one not harmonized yet is
    an error, as for Babel."""
    if icd10cm_nodes is None:
        return {}
    if not icd10cm_nodes.exists():
        raise FileNotFoundError(f"{icd10cm_nodes}: icd10cm must be harmonized before WHO's ICD-10")
    return harmonized_titles(icd10cm_nodes)


def record_of_reads(placed_by_babel: set[str], cm_titles: dict[str, str]) -> dict[str, int | str]:
    """What the output was made from, as ``READS_FILENAME`` holds it: a digest of the two facts read, so that a
    copied file is the same and a changed one is not."""
    lines = [*sorted(placed_by_babel), *(f"{code}\t{title}" for code, title in sorted(cm_titles.items()))]
    return {
        "codes_babel_places": len(placed_by_babel),
        "icd10cm_titles": len(cm_titles),
        "sha256": hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest(),
    }


def output_is_current(nodes_output: Path, babel_edges: Path | None, icd10cm_nodes: Path | None) -> bool:
    """Whether the harmonized output was made from what Babel's output and icd10cm's hold now."""
    recorded = Path(nodes_output).with_name(READS_FILENAME)
    return recorded.exists() and json.loads(recorded.read_text()) == record_of_reads(
        codes_babel_places(babel_edges), icd10cm_titles(icd10cm_nodes)
    )


class ICD10Harmonizer(BaseHarmonizer):
    """Harmonizer for the code file of WHO's ICD-10."""

    def __init__(self, *args, babel_edges: Path | None = None, icd10cm_nodes: Path | None = None, **kwargs):
        """``babel_edges``: Babel's harmonized edges, read for the codes Babel already places. ``icd10cm_nodes``:
        the harmonized nodes of the build's icd10cm source, read for the twins. So this source is harmonized after
        both (see ``orchestrator.harmonization_order``). Either is None where the build lacks that source."""
        super().__init__(*args, **kwargs)
        self.babel_edges = babel_edges
        self.icd10cm_nodes = icd10cm_nodes

    def harmonize(
        self,
        nodes_output: Path,
        edges_output: Path,
        *,
        input_file: Path | None = None,
        nodes_input: Path | None = None,
        edges_input: Path | None = None,
    ):
        if not input_file:
            raise ValueError(f"{self.source_name} requires input_file (the directory holding WHO's code file)")
        input_dir = Path(input_file)
        logging.info(f"Harmonizing {self.source_name}: {input_dir} -> {nodes_output}, {edges_output}")
        record = Path(nodes_output).with_name(READS_FILENAME)
        record.unlink(missing_ok=True)  # a run that dies half way must not leave the last run's record on its files

        entries = read_codes_file(single_file(input_dir, CODES_FILE_PATTERN))
        names = code_names(entries)
        cm_titles = icd10cm_titles(self.icd10cm_nodes)
        placed_by_babel = codes_babel_places(self.babel_edges)

        nodes, edges = [], []
        stubs: dict[str, dict] = {}

        def stub(curie: str) -> str:
            stubs.setdefault(
                curie, self.create_node(curie=curie, categories=STUB_CATEGORY, provided_by=self.source_infores)
            )
            return curie

        for code, entry in entries.items():
            if who_curie(code) in placed_by_babel:
                continue
            has_twin = code in cm_titles and title_key(cm_titles[code]) == title_key(entry.title)
            matched_on = not (has_twin or titled_by_agent(code, entry))
            node = self.create_node(
                curie=who_curie(code),
                # ICD-10-CM's category for the chapter. A code moved to its twin lends that node no category.
                categories=CATEGORY_BY_CHAPTER.get(code[0], DISEASE_CATEGORY),
                provided_by=self.source_infores,
                name=names[code] if matched_on else None,
                unmatched_name=None if matched_on else names[code],
                attributes={LEAF_ATTRIBUTE: entry.leaf},
            )
            nodes.append(node)
            parent = parent_code(code, entries)
            if parent and who_curie(parent) not in placed_by_babel:
                edges.append(
                    self.create_edge(
                        subject_id=node["id"],
                        object_id=who_curie(parent),
                        predicate=SUBCLASS_OF_PREDICATE,
                        primary_ks=self.source_infores,
                        knowledge_level=KNOWLEDGE_ASSERTION,
                        agent_type=MANUAL_AGENT,
                    )
                )
            if has_twin:
                # Our inference from the two lists, not a statement of WHO's or CDC's
                edges.append(
                    self.create_edge(
                        subject_id=node["id"],
                        object_id=stub(cm_curie(code)),
                        predicate=SAME_AS_PREDICATE,
                        primary_ks=KRAKEN_SOURCE_ID,
                        knowledge_level=KNOWLEDGE_ASSERTION,
                        agent_type=AUTOMATED_AGENT,
                        supporting_sources=[self.source_infores, ICD10CM_INFORES],
                    )
                )

        logging.info(
            f"{len(entries)} codes: {len(entries) - len(nodes)} left as Babel has them, {len(nodes)} nodes, "
            f"{sum(edge['predicate'] == SAME_AS_PREDICATE for edge in edges)} of them the ICD-10-CM code of that code"
        )
        save_to_jsonl([*nodes, *stubs.values()], nodes_output, mode="w")
        save_to_jsonl(edges, edges_output, mode="w")
        record.write_text(json.dumps(record_of_reads(placed_by_babel, cm_titles)))
        logging.info(
            f"{self.source_name} harmonization complete: {len(nodes)} nodes, {len(stubs)} stubs, {len(edges)} edges"
        )
