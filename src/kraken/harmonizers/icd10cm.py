# icd10cm.py
"""ICD-10-CM -- the diagnosis codes in US patient records -- so a coded diagnosis can reach genes and drugs.

Input is a directory (config ``input_file``) holding CDC's order file, which lists every code and header of a fiscal
year, and dated copies of Mondo's SSSOM mappings and Disease Ontology's release, the sources of the disease links:

  * every code is a node with a ``biolink:subclass_of`` edge to its parent. The tree is in the code itself: a code's
    parent is the longest prefix of it that the file also lists (``T36.0X1A`` -> ``T36.0X1`` -> ``T36.0``).
  * a Mondo exact row is a ``biolink:same_as`` edge, which entity resolution merges on, unless the row may not be
    exact (``mondo_links`` says when). Then it is a ``biolink:close_match``.
  * a code with no Mondo row of its own (``I10``) gets a ``biolink:close_match`` from every Disease Ontology term
    that lists it as an xref. Never a same_as: an xref ties a code to any term classed
    under it (RIDDLE syndrome lists ``D82.8``), and a close_match stays below entity resolution's merge threshold.

Mondo's broad rows (the code is broader than the Mondo term) are left out. Only 60 of the 78 name a code in the file
and a live term, and 38 of those are a residual code ("Other specified immunodeficiencies" for one rare disease).
"""

import csv
import logging
import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import NamedTuple

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.constants import CLOSE_MATCH_PREDICATE, KNOWLEDGE_ASSERTION, MANUAL_AGENT, SAME_AS_PREDICATE
from kraken.utils.kg_io import save_to_jsonl

ORDER_FILE_PATTERN = "icd10cm-order-*.txt"
MONDO_MAPPINGS_PATTERN = "mondo.sssom.*.tsv"
DOID_FILE_PATTERN = "doid.*.obo"

ICD10CM_PREFIX = "ICD10CM"
SHORTEST_CODE_LENGTH = 3
RESIDUAL_NAME = re.compile(r"\bother\b", re.IGNORECASE)  # "Other specified diabetes mellitus": a group, not a disease

MONDO_INFORES = "infores:mondo"
MONDO_AGENT_TYPE = "manual_validation_of_automated_agent"  # what the infores catalog records for Mondo
SKOS_EXACT_MATCH = "skos:exactMatch"
OBSOLETE_LABEL_PREFIX = "obsolete "
# Exact rows reviewed as not exact. Merged, a broader code puts its other diagnoses under the term (M33.2
# "Polymyositis" under dermatomyositis), and a narrower code gives the term its own parent (aortic valve
# insufficiency under I06 "Rheumatic aortic valve diseases"). Code -> the term its row names.
NOT_EXACT_MONDO_ROWS: dict[str, str] = {
    # the code is broader than the term
    "A59.0": "MONDO:0005993",  # Trichomonas vaginitis; the code also has trichomonal prostatitis
    "B07": "MONDO:0001209",  # common wart; the code also has plantar wart
    "C50.2": "MONDO:0002067",  # female breast cancers; each of these three codes has a male child
    "C50.3": "MONDO:0002068",
    "C50.6": "MONDO:0002069",
    "D01.5": "MONDO:0004715",  # liver carcinoma in situ; the code adds gallbladder and bile ducts
    "E56.1": "MONDO:0001244",  # vitamin K deficiency hemorrhagic disease; the code needs no bleeding
    "E71.310": "MONDO:0008723",  # very long chain; the code also has long chain
    "H31.22": "MONDO:0008982",  # central areolar; the code also has generalized and peripapillary
    "H54": "MONDO:0001941",  # blindness; the code also has low vision
    "I80.1": "MONDO:0001481",  # thrombophlebitis; the code also has phlebitis
    "J94.2": "MONDO:0006783",  # hemopneumothorax; hemothorax needs no air
    "K52.81": "MONDO:0016129",  # eosinophilic gastroenteritis; the code also has gastritis
    "M33": "MONDO:0016367",  # dermatomyositis; the code also has polymyositis
    "P29.1": "MONDO:0018684",  # idiopathic neonatal atrial flutter; the code also has tachycardia and bradycardia
    "Q23": "MONDO:0019817",  # mitral valve; the code also has aortic valve
    "Q44.6": "MONDO:0000447",  # autosomal dominant polycystic liver disease; the code is any cystic liver disease
    # the code is narrower than the term: what else the term covers is in brackets
    "A18.2": "MONDO:0005831",  # lymph node tuberculosis [intrathoracic, A15.4]
    "A66.4": "MONDO:0001066",  # late yaws [gangosa, bone and joint lesions]
    "B56": "MONDO:0000940",  # trypanosomiasis [Chagas disease]
    "D57.2": "MONDO:0011382",  # sickle cell disease [Hb-SS, D57.0, and other genotypes]
    "E71.510": "MONDO:0019609",  # Zellweger spectrum disorders [neonatal adrenoleukodystrophy, E71.511]
    "H00": "MONDO:0003382",  # eyelid disorder [every other one]
    "H46.2": "MONDO:0001431",  # toxic or nutritional optic neuropathy [toxic, H46.3]
    "H53.42": "MONDO:0004758",  # scotoma [central and other scotomas]
    "H57.01": "MONDO:0001986",  # Argyll Robertson pupil [the typical, syphilitic one]
    "H65.2": "MONDO:0005892",  # otitis media with effusion [acute serous, H65.0]
    "H68": "MONDO:0004866",  # eustachian tube disorder [other ones, H69]
    "I06.1": "MONDO:0005648",  # aortic valve insufficiency [nonrheumatic, I35.1]
    "I27.81": "MONDO:0004596",  # cor pulmonale [acute]
    "J04.1": "MONDO:0005990",  # tracheitis [chronic]
    "J30": "MONDO:0003014",  # rhinitis [chronic, J31]
    "J35.01": "MONDO:0001039",  # tonsillitis [acute, J03]
    "J35.02": "MONDO:0000261",  # adenoiditis [acute]
    "K30": "MONDO:0002268",  # dyspepsia [not functional]
    "K31.82": "MONDO:0001427",  # Dieulafoy lesion [of the intestine, K63.81]
    "M42": "MONDO:0018381",  # osteochondrosis [outside the spine]
    "M75.1": "MONDO:0007028",  # rotator cuff syndrome [impingement, M75.4]
    "M77.1": "MONDO:0001875",  # epicondylitis [medial, M77.0]
    "M81": "MONDO:0005298",  # osteoporosis [with a current fracture, M80]
    "N25": "MONDO:0001343",  # impaired renal function disease [not tubular]
    "N84": "MONDO:0005079",  # polyp [outside the female genital tract]
    "P14.0": "MONDO:0700303",  # Erb palsy [not from birth injury]
    "Q96.0": "MONDO:0019499",  # Turner syndrome [mosaics and other karyotypes, Q96.1-Q96.9]
}

DOID_INFORES = "infores:disease-ontology"
DOID_XREF_LINE = "xref: ICD10CM:"
OBO_XREF = "oboInOwl:hasDbXref"

SUBCLASS_OF_PREDICATE = "biolink:subclass_of"
BILLABLE_ATTRIBUTE = "billable"
ORIGINAL_PREDICATE_ATTRIBUTE = "original_predicate"

DISEASE_CATEGORY = "biolink:Disease"
# Every category here must sit in exactly ONE entity resolution family. ER only inherits a source's category from a
# single-family node, so an id Babel doesn't know that is typed with a bridge such as ClinicalFinding falls through
# to NamedThing, a wildcard that can merge on name with anything. Babel knows no ICD10CM id.
CATEGORY_BY_CHAPTER: dict[str, str] = {
    "R": "biolink:PhenotypicFeature",  # symptoms, signs and abnormal findings
    # Injury and poisoning: Biolink's disease excludes "a direct result of physical injury"
    "S": "biolink:DiseaseOrPhenotypicFeature",
    "T": "biolink:DiseaseOrPhenotypicFeature",
    # External causes ("Pecked by chicken, initial encounter") happen to a patient, they aren't conditions
    "V": "biolink:Event",
    "W": "biolink:Event",
    "X": "biolink:Event",
    "Y": "biolink:Event",
    "Z": "biolink:ClinicalEntity",  # "Encounter for ...", "Personal history of ...": facts about care
}


class OrderEntry(NamedTuple):
    billable: bool
    name: str


class MondoRow(NamedTuple):
    mondo_id: str
    mondo_label: str
    predicate: str
    target: str  # for ICD-10-CM, a dotted code or a range such as ICD10CM:I10-I16


class MondoLink(NamedTuple):
    mondo_id: str
    mondo_label: str
    code: str
    predicate: str


class DoidTerm(NamedTuple):
    doid_id: str
    name: str


def read_order_file(path: Path) -> dict[str, OrderEntry]:
    """``code without its dot -> entry``. Layout from CDC's icd10OrderFiles.pdf: the code at 7-13, the billable flag
    at 15 and the long name from 78."""
    entries: dict[str, OrderEntry] = {}
    with open(path, encoding="ascii") as f:
        for line in f:
            code, flag, name = line[6:13].strip(), line[14:15], line[77:].strip()
            if flag not in ("0", "1") or not code or not name:
                raise ValueError(f"{path}: not an order file line: {line!r}")
            entries[code] = OrderEntry(billable=flag == "1", name=name)
    return entries


def read_mondo_rows(path: Path) -> Iterator[MondoRow]:
    """Every row of Mondo's mappings but those of obsolete terms: an obsolete term is no disease to link a code to."""
    with open(path, encoding="utf-8", newline="") as f:
        lines = (line for line in f if not line.startswith("#"))
        for row in csv.DictReader(lines, delimiter="\t"):
            if not row["subject_label"].startswith(OBSOLETE_LABEL_PREFIX):
                yield MondoRow(row["subject_id"], row["subject_label"], row["predicate_id"], row["object_id"])


def read_doid_xrefs(path: Path) -> dict[str, set[DoidTerm]]:
    """``code without its dot -> terms`` from the ICD10CM xrefs of live Disease Ontology terms. An OBO file is
    stanzas split by blank lines; a term's xrefs are lines like ``xref: ICD10CM:I10``."""
    terms_by_code: dict[str, set[DoidTerm]] = defaultdict(set)
    for stanza in path.read_text(encoding="utf-8").split("\n\n"):
        lines = stanza.splitlines()
        if not lines or lines[0] != "[Term]" or "is_obsolete: true" in lines:
            continue
        fields = {key: value for key, _sep, value in (line.partition(": ") for line in lines) if key in ("id", "name")}
        term = DoidTerm(fields["id"], fields["name"])
        for line in lines:
            if line.startswith(DOID_XREF_LINE):
                terms_by_code[line.removeprefix(DOID_XREF_LINE).replace(".", "")].add(term)
    return terms_by_code


def dotted(code: str) -> str:
    """The code as it is written everywhere but CDC's files: a dot after the third character."""
    return code if len(code) <= SHORTEST_CODE_LENGTH else f"{code[:SHORTEST_CODE_LENGTH]}.{code[SHORTEST_CODE_LENGTH:]}"


def cm_curie(code: str) -> str:
    return f"{ICD10CM_PREFIX}:{dotted(code)}"


def parent_code(code: str, codes: dict[str, OrderEntry]) -> str | None:
    """The longest proper prefix of ``code`` the file lists. A 3-character code is a root."""
    for end in range(len(code) - 1, SHORTEST_CODE_LENGTH - 1, -1):
        if code[:end] in codes:
            return code[:end]
    if len(code) > SHORTEST_CODE_LENGTH:
        raise ValueError(f"ICD-10-CM code {code} has no parent in the order file")
    return None


def _name_key(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", text.lower())


def mondo_links(rows: Iterable[MondoRow], codes: dict[str, OrderEntry]) -> list[MondoLink]:
    """A same_as only where Mondo's exact row can be merged on. Otherwise a close_match: the disease stays linked, the
    code keeps its name, and a name match decides whether it merges. Not every exact row is exact. Not merged on:

      * a term with several targets, since they can't all be it. 13 of the 15 terms with two codes hold a parent and its
        child (``K28`` and ``K28.0``, "Acute gastrojejunal ulcer with hemorrhage"), and nothing in the row tells which
        is the term: it is the parent for gastrojejunal ulcer, the child ``F53.1`` for postpartum psychosis. A range
        counts as a target though it is never linked: ``I10-I16`` is hypertensive disorder, the same term's ``I15`` is
        not.
      * a residual code: ``E13`` "Other specified diabetes mellitus" is exact to cystic fibrosis-related diabetes.
      * a code named exactly like another term: ``D70`` "Neutropenia" is exact to agranulocytosis, though Mondo has
        neutropenia. Its name takes it there.
      * a row reviewed as not exact (``NOT_EXACT_MONDO_ROWS``).
    """
    label_keys: set[str] = set()
    rows_by_term: dict[str, list[MondoRow]] = defaultdict(list)
    for row in rows:
        label_keys.add(_name_key(row.mondo_label))
        if row.predicate == SKOS_EXACT_MATCH and row.target.startswith(f"{ICD10CM_PREFIX}:"):
            rows_by_term[row.mondo_id].append(row)

    links = []
    for mondo_id, term_rows in rows_by_term.items():
        for row in term_rows:
            dotted_code = row.target.split(":", 1)[1]
            code = dotted_code.replace(".", "")
            if code not in codes:  # a range, or a code this fiscal year lacks
                continue
            name = codes[code].name
            named_like_another_term = _name_key(name) in label_keys and _name_key(name) != _name_key(row.mondo_label)
            reviewed_not_exact = NOT_EXACT_MONDO_ROWS.get(dotted_code) == mondo_id
            mergeable = (
                len(term_rows) == 1
                and not RESIDUAL_NAME.search(name)
                and not named_like_another_term
                and not reviewed_not_exact
            )
            predicate = SAME_AS_PREDICATE if mergeable else CLOSE_MATCH_PREDICATE
            links.append(MondoLink(mondo_id, row.mondo_label, code, predicate))
    held_back = [dotted(link.code) for link in links if link.predicate == CLOSE_MATCH_PREDICATE]
    logging.info(f"{len(held_back)} of {len(links)} usable Mondo exact rows are close_match, not same_as: {held_back}")
    if fixed := sorted(NOT_EXACT_MONDO_ROWS.items() - {(dotted(link.code), link.mondo_id) for link in links}):
        logging.warning(f"Mondo no longer has these exact rows, so drop them from NOT_EXACT_MONDO_ROWS: {fixed}")
    return links


class ICD10CMHarmonizer(BaseHarmonizer):
    """Harmonizer for CDC's ICD-10-CM order file, and Mondo's and Disease Ontology's links to it."""

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
            raise ValueError(f"{self.source_name} requires input_file (the directory holding CDC's and Mondo's files)")
        input_dir = Path(input_file)
        logging.info(f"Harmonizing {self.source_name}: {input_dir} -> {nodes_output}, {edges_output}")

        codes = read_order_file(self._single_file(input_dir, ORDER_FILE_PATTERN))
        links = mondo_links(read_mondo_rows(self._single_file(input_dir, MONDO_MAPPINGS_PATTERN)), codes)
        merging_codes = {link.code for link in links if link.predicate == SAME_AS_PREDICATE}

        nodes: dict[str, dict] = {}
        edges: list[dict] = []
        for code, entry in codes.items():
            node = self._code_node(code, entry, merges_through_mondo=code in merging_codes)
            nodes[node["id"]] = node
            parent = parent_code(code, codes)
            if parent:
                edges.append(
                    self.create_edge(
                        subject_id=node["id"],
                        object_id=cm_curie(parent),
                        predicate=SUBCLASS_OF_PREDICATE,
                        primary_ks=self.source_infores,
                        knowledge_level=KNOWLEDGE_ASSERTION,
                        agent_type=MANUAL_AGENT,
                    )
                )

        for link in links:
            mondo_node = self.create_node(
                curie=link.mondo_id, categories=DISEASE_CATEGORY, provided_by=self.source_infores, name=link.mondo_label
            )
            nodes.setdefault(mondo_node["id"], mondo_node)
            edges.append(
                self.create_edge(
                    subject_id=link.mondo_id,
                    object_id=cm_curie(link.code),
                    predicate=link.predicate,
                    primary_ks=MONDO_INFORES,
                    knowledge_level=KNOWLEDGE_ASSERTION,
                    agent_type=MONDO_AGENT_TYPE,
                    attributes={ORIGINAL_PREDICATE_ATTRIBUTE: SKOS_EXACT_MATCH},
                )
            )

        # A code's own Mondo row wins: Disease Ontology's term is then nearly always in the same cluster. A parent's
        # row doesn't: E74.05 "LAMP2 deficiency" keeps Danon disease, though E74.0 reaches "disorder of glycogen
        # metabolism".
        mondo_codes = {link.code for link in links}
        doid_terms = read_doid_xrefs(self._single_file(input_dir, DOID_FILE_PATTERN))
        for code, terms in sorted(doid_terms.items()):
            if code not in codes or code in mondo_codes or RESIDUAL_NAME.search(codes[code].name):
                continue  # a range, a code this fiscal year lacks, or a residual group
            for term in sorted(terms):
                # Named, like Mondo's stubs: a term newer than our Babel would otherwise be a nameless node
                doid_node = self.create_node(
                    curie=term.doid_id, categories=DISEASE_CATEGORY, provided_by=self.source_infores, name=term.name
                )
                nodes.setdefault(doid_node["id"], doid_node)
                edges.append(
                    self.create_edge(
                        subject_id=term.doid_id,
                        object_id=cm_curie(code),
                        predicate=CLOSE_MATCH_PREDICATE,
                        primary_ks=DOID_INFORES,
                        knowledge_level=KNOWLEDGE_ASSERTION,
                        agent_type=MANUAL_AGENT,
                        attributes={ORIGINAL_PREDICATE_ATTRIBUTE: OBO_XREF},
                    )
                )

        logging.info(f"Saving {len(nodes)} nodes and {len(edges)} edges")
        save_to_jsonl(nodes.values(), nodes_output, mode="w")
        save_to_jsonl(edges, edges_output, mode="w")
        logging.info(f"{self.source_name} harmonization complete: {len(nodes)} nodes, {len(edges)} edges")

    @staticmethod
    def _single_file(input_dir: Path, pattern: str) -> Path:
        matches = sorted(input_dir.glob(pattern))
        if len(matches) != 1:
            raise ValueError(f"expected one {pattern} in {input_dir}, found {[m.name for m in matches]}")
        return matches[0]

    def _code_node(self, code: str, entry: OrderEntry, *, merges_through_mondo: bool) -> dict:
        """A code that merges through Mondo carries its name as a synonym only. Entity resolution merges on an equal
        name, so a name here would join whatever else has it to the Mondo cluster through this code, and a wrong row
        would take more than the code: ``J94.2`` "Hemothorax", exact to hemopneumothorax, would have taken the
        hemothorax phenotype HP:0012151 with it. Mondo names the merged node anyway. And that code is a disease
        whatever its chapter says, or the family guardrail would cut the same_as (``Z15``)."""
        return self.create_node(
            curie=cm_curie(code),
            categories=DISEASE_CATEGORY if merges_through_mondo else CATEGORY_BY_CHAPTER.get(code[0], DISEASE_CATEGORY),
            provided_by=self.source_infores,
            name=None if merges_through_mondo else entry.name,
            synonyms={entry.name},
            attributes={BILLABLE_ATTRIBUTE: entry.billable},
        )
