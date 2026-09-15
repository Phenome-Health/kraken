# babel.py
"""Babel -- the identifier-equivalence build behind the SRI Node Normalizer -- ingested as a source.

Babel publishes, per Biolink type, a "compendium" of cliques: sets of identifiers it considers the same entity,
each identifier with its OWN label, descriptions and taxa. KRAKEN used to get this by querying the Node
Normalizer API during entity resolution; reading the release files instead makes it versioned, reproducible and
free of API calls, and lets us choose the conflations ourselves.

REPRESENTATION. Every identifier is its own node (id, name, description, taxon, the clique's category). A
clique is NOT folded into equivalent_ids -- that would keep only one name per clique, and entity resolution
needs each id's own name. Instead each clique is a star of ``biolink:same_as`` edges from its preferred id
(Babel lists it first) to every other member; entity resolution regroups the star back into a clique (see
entity_resolution.build). Babel's two conflations are handled differently, on purpose:

  * GENE/PROTEIN conflation is ON: ``gene biolink:same_as protein`` edges, marked as the conflation in their
    BABEL_RELATION_ATTRIBUTE, which entity resolution treats as merge evidence (the Node Normalizer's ``conflate=true``
    did the same). Where entity resolution keeps a gene and its protein apart anyway, integration turns the edge
    into a close_match, like any other Babel same_as between clusters.
  * DRUG/CHEMICAL conflation is OFF: a drug product, its salt forms and its active compound stay separate
    entities. What Babel knows about how they relate is kept as typed edges instead (RELATION_PREDICATES), plus
    a ``biolink:close_match`` from each conflation group's preferred id to any member no typed relation reaches.
    None of these are merge evidence.

SCOPE. The full release is ~414M identifiers, ten times KRAKEN. Nearly all of the excess is three things, each
filtered by a rule that looks only at the clique itself -- never at what other sources contain, so the result
doesn't depend on build order:
  * Gene / Protein cliques from organisms of no biomedical interest: kept only if a clique taxon is in
    ncbigene_taxon_allowlist.py (after rolling up to species -- the same list NCBI Gene is scoped by).
  * SmallMolecule / MolecularMixture cliques made only of structure-registry ids (PubChem, InChIKey, CAS, ChEMBL):
    ~117M structures no curated vocabulary names -- vendor compounds, screening libraries, patent examples. Kept
    only if some other identifier (HMDB, ChEBI, UNII, MeSH, DrugBank, ...) is present: 1.9M of SmallMolecule's
    232M ids. Counting CAS or ChEMBL as curation instead would keep 10M or 21M.
  * Publication: excluded entirely (PMIDs aren't entities we resolve).
And in every compendium, a clique that is ONE identifier with no label and no taxon is dropped: it says nothing
beyond "this id exists" (2.5M of them are nameless Ensembl genes of unknown organism). An id like that which
another source actually uses still becomes a node through that source.
"""

import csv
import json
import logging
import re
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

import jsonlines

from kraken.harmonizers.base import BaseHarmonizer
from kraken.harmonizers.helpers.ncbigene_taxon_allowlist import TAXON_ALLOWLIST
from kraken.utils.constants import (
    AUTOMATED_AGENT,
    BABEL_RELATION_ATTRIBUTE,
    CLOSE_MATCH_PREDICATE,
    DRUG_CHEMICAL_CONFLATION_RELATION,
    GENE_PROTEIN_CONFLATION_RELATION,
    KNOWLEDGE_ASSERTION,
    SAME_AS_PREDICATE,
)
from kraken.utils.taxonomy import TaxonNormalizer

# --- Input layout (the release's own directory names, plus the two files we add alongside) ---
COMPENDIA_DIR = "compendia"
CONFLATION_DIR = "conflation"
COMPENDIUM_SUFFIX = ".txt"
GENE_PROTEIN_CONFLATION_FILENAME = "GeneProtein.txt"
DRUG_CHEMICAL_CONFLATION_FILENAME = "DrugChemical.txt"
# Babel ships its typed drug/chemical relations only inside duckdb/Concord.parquet (4.6 GB). The release
# directory's extract_drugchemical_concords.sh pulls just those rows into this TSV (subj, pred, obj).
DRUG_CHEMICAL_RELATIONS_FILENAME = "DrugChemical_concords.tsv"
# NCBI's taxonomy, for rolling taxa up to species before the allowlist check (see NCBIGeneHarmonizer).
TAXDUMP_FILENAME = "taxdump.tar.gz"

# --- Scope (see the module docstring) ---
EXCLUDED_COMPENDIA = frozenset({"Publication"})
TAXON_FILTERED_COMPENDIA = frozenset({"Gene", "Protein"})
STRUCTURE_FILTERED_COMPENDIA = frozenset({"SmallMolecule", "MolecularMixture"})
# Identifiers assigned to structures wholesale -- every deposited compound (PubChem, InChIKey), every registered
# substance (CAS), every compound in a medicinal-chemistry paper or assay (ChEMBL) -- so on their own they say
# nothing about whether anyone curates the entity.
STRUCTURE_ONLY_PREFIXES = frozenset({"PUBCHEM.COMPOUND", "INCHIKEY", "CAS", "CHEMBL.COMPOUND"})

# --- Edges ---
# Babel's drug/chemical relations (RxNorm's, via UMLS) read subject -> object: "doxepin 100 MG Oral Capsule
# has_active_ingredient doxepin hydrochloride". Biolink 4.2.5 has no plain has_ingredient/has_constituent, so an
# ingredient that isn't known to be the ACTIVE one is has_part. The rest say two ids name nearly the same thing
# (a brand and its generic, two dose forms of one drug, an RxNorm ingredient and its PubChem compound).
RELATION_PREDICATES: dict[str, str] = {
    "has_active_ingredient": "biolink:has_active_ingredient",
    "has_precise_active_ingredient": "biolink:has_active_ingredient",
    "has_ingredient": "biolink:has_part",
    "has_precise_ingredient": "biolink:has_part",
    "consists_of": "biolink:has_part",
    "tradename_of": CLOSE_MATCH_PREDICATE,
    "has_form": CLOSE_MATCH_PREDICATE,
    "linked": CLOSE_MATCH_PREDICATE,
}
INFORMATION_CONTENT_ATTRIBUTE = "information_content"

# Clique lines end with '"taxa": [...]}' -- read straight off the line so the ~95% of Gene/Protein cliques
# outside the allowlist are rejected without parsing their JSON.
TAXA_TAIL_PATTERN = re.compile(r'"taxa": \[([^\]]*)\]\}\s*$')
IDENTIFIER_PREFIX_PATTERN = re.compile(r'"i": "([^":]+):')


class BabelHarmonizer(BaseHarmonizer):
    """Harmonizer for a Babel release directory (https://stars.renci.org/var/babel_outputs/<release>/).

    Input is a DIRECTORY (config `input_file`) holding `compendia/*.txt`, `conflation/GeneProtein.txt`,
    `conflation/DrugChemical.txt`, `conflation/DrugChemical_concords.tsv` and `taxdump.tar.gz`.
    """

    # Keep a Gene/Protein clique that names no taxon at all. The allowlist can't judge it, and dropping it
    # would lose e.g. a gene concept a curated vocabulary names without an organism.
    keep_untaxoned_gene_protein_cliques: bool = True

    def __init__(self, biolink_client, source_id: str, **kwargs):
        super().__init__(biolink_client, source_id, **kwargs)
        self.taxonomy: TaxonNormalizer | None = None
        self._taxon_allowed: dict[str, bool] = {}
        self._species_curie: dict[str, str] = {}
        # Preferred ids of kept Gene / Protein cliques -- the conflation file lists only preferred ids.
        self.gene_leaders: set[str] = set()
        self.protein_leaders: set[str] = set()
        # Drug/chemical relation endpoints (a small set) -> the preferred id of the kept clique holding them.
        self.drug_chemical_ids: set[str] = set()
        self.drug_chemical_leader: dict[str, str] = {}
        self.stats: dict[str, Counter] = defaultdict(Counter)

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
            raise ValueError(f"{self.source_name} requires input_file (the Babel release directory)")
        release_dir = Path(input_file)
        compendia = self._find_compendia(release_dir)
        conflation_dir = release_dir / CONFLATION_DIR
        for required in (
            conflation_dir / GENE_PROTEIN_CONFLATION_FILENAME,
            conflation_dir / DRUG_CHEMICAL_CONFLATION_FILENAME,
            conflation_dir / DRUG_CHEMICAL_RELATIONS_FILENAME,
            release_dir / TAXDUMP_FILENAME,
        ):
            if not required.is_file():
                raise FileNotFoundError(f"{self.source_name}: missing {required}")

        logging.info(f"Harmonizing {self.source_name}: {len(compendia)} compendia in {release_dir}")
        self.taxonomy = TaxonNormalizer(release_dir / TAXDUMP_FILENAME, include_names=False)
        allowed_species = frozenset(self.taxonomy.to_species(tax_id) for tax_id in TAXON_ALLOWLIST)
        relations = self._load_relations(conflation_dir / DRUG_CHEMICAL_RELATIONS_FILENAME)
        groups = self._load_groups(conflation_dir / DRUG_CHEMICAL_CONFLATION_FILENAME)
        self.drug_chemical_ids = {curie for row in relations for curie in (row[0], row[2])}
        self.drug_chemical_ids.update(curie for group in groups for curie in group)

        with jsonlines.open(nodes_output, "w") as nodes, jsonlines.open(edges_output, "w") as edges:
            for compendium_path in compendia:
                self._harmonize_compendium(compendium_path, allowed_species, nodes, edges)
            self._write_gene_protein_conflations(conflation_dir / GENE_PROTEIN_CONFLATION_FILENAME, edges)
            self._write_drug_chemical_edges(relations, groups, edges)

        self._log_run_summary()

    # ------------------------------- inputs -------------------------------

    def _find_compendia(self, release_dir: Path) -> list[Path]:
        compendia_dir = release_dir / COMPENDIA_DIR
        if not compendia_dir.is_dir():
            raise FileNotFoundError(f"{self.source_name}: no {COMPENDIA_DIR}/ directory in {release_dir}")
        compendia = sorted(
            path
            for path in compendia_dir.iterdir()
            if path.suffix == COMPENDIUM_SUFFIX and path.stem not in EXCLUDED_COMPENDIA
        )
        if not compendia:
            raise FileNotFoundError(f"{self.source_name}: no compendium files in {compendia_dir}")
        return compendia

    @staticmethod
    def _load_relations(path: Path) -> list[tuple[str, str, str]]:
        with open(path, newline="") as file:
            reader = csv.reader(file, delimiter="\t")
            next(reader, None)  # header: subj, pred, obj
            return [(row[0], row[1], row[2]) for row in reader if len(row) == 3]

    @staticmethod
    def _load_groups(path: Path) -> list[list[str]]:
        with open(path) as file:
            return [json.loads(line) for line in file if line.strip()]

    # ------------------------------- compendia -------------------------------

    def _harmonize_compendium(self, path: Path, allowed_species: frozenset[str], nodes, edges) -> None:
        compendium = path.stem
        stats = self.stats[compendium]
        logging.info(f"Streaming {path.name}..")
        with open(path) as file:
            for line_num, line in enumerate(file, 1):
                if line_num % 10_000_000 == 0:
                    logging.info(f"    at {line_num} cliques of {path.name} ({stats['cliques_kept']} kept)")
                stats["cliques_read"] += 1
                if compendium in TAXON_FILTERED_COMPENDIA and not self._taxa_allowed(line, allowed_species):
                    stats["cliques_dropped_taxon"] += 1
                    continue
                if compendium in STRUCTURE_FILTERED_COMPENDIA and self._structure_only(line):
                    stats["cliques_dropped_structure_only"] += 1
                    continue
                clique = json.loads(line)
                if self._is_empty_singleton(clique):
                    stats["cliques_dropped_empty_singleton"] += 1
                    continue
                members = self._write_clique(clique, nodes, edges)
                if not members:
                    continue
                stats["cliques_kept"] += 1
                stats["ids_kept"] += len(members)
                if compendium == "Gene":
                    self.gene_leaders.add(members[0])
                elif compendium == "Protein":
                    self.protein_leaders.add(members[0])

    def _taxa_allowed(self, line: str, allowed_species: frozenset[str]) -> bool:
        match = TAXA_TAIL_PATTERN.search(line)
        if match is None:  # not the expected layout -- fall back to parsing the clique
            taxa = json.loads(line).get("taxa") or []
        else:
            taxa = [taxon.strip().strip('"') for taxon in match.group(1).split(",") if taxon.strip()]
        if not taxa:
            return self.keep_untaxoned_gene_protein_cliques
        for taxon in taxa:
            allowed = self._taxon_allowed.get(taxon)
            if allowed is None:
                allowed = self.taxonomy.to_species(taxon.removeprefix("NCBITaxon:")) in allowed_species
                self._taxon_allowed[taxon] = allowed
            if allowed:
                return True
        return False

    @staticmethod
    def _is_empty_singleton(clique: dict[str, Any]) -> bool:
        identifiers = clique.get("identifiers") or []
        return len(identifiers) == 1 and not identifiers[0].get("l") and not identifiers[0].get("t")

    @staticmethod
    def _structure_only(line: str) -> bool:
        return all(prefix in STRUCTURE_ONLY_PREFIXES for prefix in IDENTIFIER_PREFIX_PATTERN.findall(line))

    def _write_clique(self, clique: dict[str, Any], nodes, edges) -> list[str]:
        """Write one node per identifier and a same_as star from the preferred (first) id. Returns the member
        ids as written (normalized), preferred id first."""
        category = clique.get("type")
        if not category:
            self.stats["_all"]["cliques_without_type"] += 1
            return []
        members: list[str] = []
        for position, identifier in enumerate(clique.get("identifiers") or []):
            raw_curie = identifier.get("i")
            if not raw_curie:
                continue
            attributes = {}
            if position == 0 and clique.get("ic") is not None:
                attributes[INFORMATION_CONTENT_ATTRIBUTE] = clique["ic"]
            node = self.create_node(
                curie=raw_curie,
                categories=category,
                provided_by=self.source_infores,
                name=identifier.get("l") or None,
                description=next((d for d in identifier.get("d") or [] if d), None),
                taxon=[self._species(taxon) for taxon in identifier.get("t") or []],
                attributes=attributes,
            )
            nodes.write(node)
            curie = node["id"]
            members.append(curie)
            if raw_curie in self.drug_chemical_ids:
                self.drug_chemical_leader[raw_curie] = members[0]
        hub = members[0] if members else None
        for member in members[1:]:
            if member != hub:
                edges.write(self._edge(hub, SAME_AS_PREDICATE, member))
        return members

    def _species(self, taxon: str) -> str:
        """A taxon CURIE rolled up to species rank (cached; there are far fewer taxa than ids)."""
        species = self._species_curie.get(taxon)
        if species is None:
            species = self.taxonomy.to_species(taxon)
            self._species_curie[taxon] = species
        return species

    # ------------------------------- conflations -------------------------------

    def _write_gene_protein_conflations(self, path: Path, edges) -> None:
        """``gene same_as protein`` for each GeneProtein conflation whose cliques we kept."""
        stats = self.stats["GeneProtein"]
        with open(path) as file:
            for line in file:
                if not line.strip():
                    continue
                ids = [self.normalize_curie(curie) for curie in json.loads(line)]
                genes = [curie for curie in ids if curie in self.gene_leaders]
                proteins = [curie for curie in ids if curie in self.protein_leaders]
                if not (genes and proteins):
                    stats["groups_without_kept_gene_and_protein"] += 1
                    continue
                for gene in genes:
                    for protein in proteins:
                        edges.write(
                            self._edge(
                                gene,
                                SAME_AS_PREDICATE,
                                protein,
                                {BABEL_RELATION_ATTRIBUTE: GENE_PROTEIN_CONFLATION_RELATION},
                            )
                        )
                        stats["gene_protein_same_as_edges"] += 1

    def _write_drug_chemical_edges(self, relations: list[tuple[str, str, str]], groups: list[list[str]], edges):
        """Typed edges for Babel's drug/chemical relations, then a close_match from each conflation group's
        preferred id to every member those relations don't connect it to."""
        stats = self.stats["DrugChemical"]
        leader = self.drug_chemical_leader
        adjacency: dict[str, set[str]] = defaultdict(set)  # between clique preferred ids, undirected
        for subject, relation, object_ in relations:
            predicate = RELATION_PREDICATES.get(relation)
            if predicate is None:
                stats[f"relations_unmapped:{relation}"] += 1
                continue
            if subject not in leader or object_ not in leader:
                stats["relations_endpoint_not_kept"] += 1
                continue
            subject_leader, object_leader = leader[subject], leader[object_]
            if subject_leader == object_leader:
                stats["relations_within_one_clique"] += 1
                continue
            adjacency[subject_leader].add(object_leader)
            adjacency[object_leader].add(subject_leader)
            edges.write(
                self._edge(
                    self.normalize_curie(subject),
                    predicate,
                    self.normalize_curie(object_),
                    {BABEL_RELATION_ATTRIBUTE: relation},
                )
            )
            stats[f"edges:{relation}"] += 1

        for group in groups:
            kept = list(dict.fromkeys(leader[curie] for curie in group if curie in leader))
            if len(kept) < 2:
                continue
            head, members = kept[0], set(kept[1:])
            reached = self._reachable(head, adjacency, set(kept))
            for member in kept[1:]:
                if member not in reached:
                    edges.write(
                        self._edge(
                            head,
                            CLOSE_MATCH_PREDICATE,
                            member,
                            {BABEL_RELATION_ATTRIBUTE: DRUG_CHEMICAL_CONFLATION_RELATION},
                        )
                    )
                    stats["conflation_close_match_edges"] += 1
            stats["conflation_members_reached_by_relations"] += len(members & reached)

    @staticmethod
    def _reachable(start: str, adjacency: dict[str, set[str]], within: set[str]) -> set[str]:
        """Ids reachable from ``start`` through relation edges, staying inside one conflation group."""
        seen, queue = {start}, deque([start])
        while queue:
            for neighbour in adjacency.get(queue.popleft(), ()):
                if neighbour in within and neighbour not in seen:
                    seen.add(neighbour)
                    queue.append(neighbour)
        return seen

    # ------------------------------- helpers -------------------------------

    def _edge(self, subject: str, predicate: str, object_: str, attributes: dict | None = None) -> dict:
        return self.create_edge(
            subject_id=subject,
            object_id=object_,
            predicate=predicate,
            primary_ks=self.source_infores,
            knowledge_level=KNOWLEDGE_ASSERTION,
            agent_type=AUTOMATED_AGENT,
            attributes=attributes,
        )

    def normalize_curie(self, curie: str) -> str:
        """As the base class, but without remembering ids biomapper2 leaves unchanged -- which for Babel is
        nearly all of them. The base cache would otherwise hold every one of Babel's tens of millions of ids."""
        if curie in self.normalized_id_map:
            return self.normalized_id_map[curie]
        prefix = curie.split(":", 1)[0]
        failures_before = self._unnormalized_count(prefix)
        normalized = super().normalize_curie(curie)
        # Forget only a SUCCESSFUL no-op. An id biomapper2 couldn't handle stays cached, so it is tallied once in
        # the normalization report rather than once per occurrence.
        if normalized == curie and self._unnormalized_count(prefix) == failures_before:
            self.normalized_id_map.pop(curie, None)
        return normalized

    def _unnormalized_count(self, prefix: str) -> int:
        return sum(
            tally[prefix]["count"]
            for tally in (self.unrecognized_vocab_prefixes, self.invalid_id_prefixes)
            if prefix in tally
        )

    def _log_run_summary(self) -> None:
        total_ids = sum(stats["ids_kept"] for stats in self.stats.values())
        total_cliques = sum(stats["cliques_kept"] for stats in self.stats.values())
        logging.info(f"{self.source_name}: kept {total_ids} identifiers in {total_cliques} cliques")
        for name, stats in sorted(self.stats.items()):
            if stats:
                logging.info(f"    {name}: {dict(stats)}")
