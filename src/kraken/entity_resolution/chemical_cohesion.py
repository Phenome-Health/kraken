"""Put back the chemical Babel cliques that clustering split (compare ``gene_protein_cohesion``).

Only a clique's own members move. Every other id stays in its cluster: a clique says which ids are one chemical,
not where the rest of their cluster belongs (joining whole clusters took lactic acid's ATC code to RefMet's
"Lactic acid").

A clique holding an InChIKey is one structure, but clustering can leave it in pieces (see ``clustering``: a clique
starts whole only within the size cap, and an id Babel has in two cliques starts alone). Under label propagation 12
of cortisone's 16 ids settled with the InChIKey of cortisone with no stereo, while its own InChIKey stayed with the
RefMet and LIPID MAPS lists that name it. No guardrail sees it, because each piece holds one InChIKey. When most of
a clique sits with an InChIKey of the same skeleton, another stereo form, those members join their own. When only a
few do, evidence of their own holds them:
GtoPdb's "L-tryptophan" sits in Babel's no-stereo clique, and its name puts it with L-tryptophan. Without a second
InChIKey there is no conflict between structures to settle, and nothing checks a member's stereo: ChEBI's
"(+)-alloaromadendrene" is in the clique of the (-) form's InChIKey. And another skeleton is another compound: MeSH's
GM3 ganglioside, a class, sits with one GM3 and is in the clique of another. Nor do members move when what stays of
their cluster may belong with them: when it would fall into pieces or leave one id alone, cut from the members that
tied it there (ibuprofen's SMILES is exactly its InChIKey), or when it keeps a SMILES, a structure no guardrail weighs
against the InChIKeys (LIPID MAPS lists all-trans-neoxanthin's InChIKey with a SMILES of the stereoisomer whose clique
would leave).

Babel also puts a UMLS concept of type T116 "Amino Acid, Peptide, or Protein" in its Protein compendium, with its MeSH
id. When that MeSH id is also in a SmallMolecule clique, the one-branch guardrail keeps the UMLS concept out of the
chemical it names (N-acetyl-S-allylcysteine). Such a clique joins the MeSH id's cluster, typed as the chemical, only
when every Protein-typed member is a UMLS, NCIT or MeSH concept that no gene or protein database id holds as a protein,
sits alone or with its clique only, and has the exact name of a member of the MeSH id's SmallMolecule clique other than
PubChem's. The shared id alone is not enough: a DrugBank link put methotrexate's brands with methopterin, MeSH's "SC002"
put an antibody-drug conjugate with a cephalosporin, and "Fibroblast Growth Factor" with one peptide's structure. Nor is
PubChem's name, which can be the MeSH name linked to the record whatever its structure ("Thioredoxins" on ergothioneine,
"Erabutoxin C", a snake toxin, on triethylphosphine sulfide), so it only repeats that link. And an id beside the
concept is there by evidence tying it to the concept, not to the chemical (LOINC's "Ristocetin", by its UMLS
equivalence).
"""

from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Callable, Collection, Iterable, Mapping
from dataclasses import dataclass, field, replace
from pathlib import Path

from kraken.entity_resolution.gene_protein_cohesion import read_cliques
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo, cluster_violations
from kraken.entity_resolution.name_sim import name_keys
from kraken.utils.constants import INCHIKEY_PREFIX, SMILES_PREFIX

PROTEIN_CATEGORY = "biolink:Protein"
GENE_CATEGORY = "biolink:Gene"
# Babel gives these two types only to an id with a structure; ChemicalEntity also holds MeSH's "POTE-2alpha protein".
STRUCTURE_CATEGORIES = frozenset({"biolink:SmallMolecule", "biolink:MolecularMixture"})
PROTEIN_DATABASE_PREFIXES = frozenset({"UniProtKB", "PR", "ENSEMBL", "CHEMBL.TARGET"})
# Where a Protein type can be UMLS's T116 on a chemical, and the one link between the compendia that held up.
RETYPABLE_PREFIXES = frozenset({"UMLS", "NCIT", "MESH"})
LINK_PREFIX = "MESH"
PUBCHEM_PREFIX = "PUBCHEM.COMPOUND"


@dataclass
class CliqueRepairs:
    replaced: list[list[str]] = field(default_factory=list)  # clusters that stop existing
    remade: list[list[str]] = field(default_factory=list)  # what their ids form instead, each guardrail-valid
    retyped: dict[str, tuple[str, ...]] = field(default_factory=dict)  # protein-typed id -> its chemical's categories


def _prefix(curie: str) -> str:
    return curie.split(":", 1)[0]


def _is_inchikey(curie: str) -> bool:
    return curie.startswith(f"{INCHIKEY_PREFIX}:")


def _is_smiles(curie: str) -> bool:
    return curie.startswith(f"{SMILES_PREFIX}:")


def _connectivity(inchikey: str) -> str:
    """An InChIKey's first block: the skeleton, without stereo, isotopes or charge."""
    return inchikey.split(":", 1)[1].split("-", 1)[0]


def plan_chemical_clique_repairs(
    curie_to_cluster: Mapping[str, int],
    cliques_path: Path,
    categories: Mapping[str, tuple[str, ...]],
    info_of: Callable[[str], NodeInfo],
    label_of: Callable[[str], str | None],
    pairs_among: Callable[[set[str]], Iterable[tuple[str, str]]],
    guardrail_config: GuardrailConfig,
    clique_cap: int,
    outliers: Collection[str] = (),
) -> CliqueRepairs:
    """What to change, without changing anything: the caller owns the clusters and their tallies. ``label_of`` is an
    id's own name, and ``pairs_among`` the match-graph pairs between the given ids."""

    def gene_or_protein_record(curie: str) -> bool:
        return _prefix(curie) in PROTEIN_DATABASE_PREFIXES or GENE_CATEGORY in categories.get(curie, ())

    def typed_protein(curie: str) -> bool:
        # An id in two Babel compendia has the type of the one read last (compendia are read in file-name order and the
        # id facts store keeps the last write): stable from build to build, but not a statement about the id.
        return categories.get(curie) == (PROTEIN_CATEGORY,)

    structure_cliques: list[list[str]] = []
    protein_cliques: list[tuple[list[str], list[str]]] = []
    in_gene_protein_clique: set[str] = set()  # only the retypable ids, a few tens of thousands
    for members in read_cliques(cliques_path):
        if any(gene_or_protein_record(m) for m in members):
            in_gene_protein_clique.update(m for m in members if _prefix(m) in RETYPABLE_PREFIXES)
            continue
        if len(members) > clique_cap:  # a star, not a clique: its members are linked only through the hub
            continue
        clustered = [m for m in members if m in curie_to_cluster and m not in outliers]
        if len({curie_to_cluster[m] for m in clustered}) < 2:
            continue
        if any(_is_inchikey(m) for m in clustered):
            structure_cliques.append(clustered)
        proteins = [m for m in clustered if typed_protein(m)]
        links = [m for m in clustered if STRUCTURE_CATEGORIES.intersection(categories.get(m, ()))]
        if proteins and links and all(_prefix(m) == LINK_PREFIX for m in links):
            protein_cliques.append((proteins, links))

    # The names of each link's SmallMolecule cliques, the chemical a protein-typed member must name exactly. Not the
    # links' own, nor PubChem's: either may be the very MeSH name in question.
    wanted_links = {m for _proteins, links in protein_cliques for m in links}
    chemical_names: dict[str, set[str]] = defaultdict(set)
    if wanted_links:
        for members in read_cliques(cliques_path):
            held = wanted_links.intersection(members)
            if held and not any(map(typed_protein, members)):
                namers = [m for m in members if m not in held and _prefix(m) != PUBCHEM_PREFIX]
                keys = {key for m in namers for key in name_keys(label_of(m), info_of(m).branches)}
                for link in held:
                    chemical_names[link] |= keys

    touched = {curie_to_cluster[m] for clique in structure_cliques for m in clique}
    touched.update(curie_to_cluster[m] for proteins, links in protein_cliques for m in proteins + links)
    members_of: dict[int, set[str]] = defaultdict(set)
    if touched:
        for curie, cluster in curie_to_cluster.items():
            if cluster in touched:
                members_of[cluster].add(curie)
    before = {cluster: sorted(members) for cluster, members in members_of.items()}
    neighbours: dict[str, set[str]] = defaultdict(set)
    for a, b in pairs_among({curie for members in members_of.values() for curie in members}):
        neighbours[a].add(b)
        neighbours[b].add(a)
    moved_to: dict[str, int] = {}
    retyped_branches: dict[str, frozenset[str]] = {}
    repairs = CliqueRepairs()

    def cluster_of(curie: str) -> int:
        return moved_to.get(curie, curie_to_cluster[curie])

    def info(curie: str) -> NodeInfo:
        node = info_of(curie)
        return replace(node, branches=retyped_branches[curie]) if curie in retyped_branches else node

    def violates(curies: Iterable[str], retype: Mapping[str, frozenset[str]] | None = None) -> bool:
        ids = sorted(curies)
        infos = {curie: info(curie) for curie in ids}
        infos.update({curie: replace(infos[curie], branches=branches) for curie, branches in (retype or {}).items()})
        return bool(cluster_violations(ids, infos, guardrail_config))

    def what_stays(movers: list[str]) -> list[set[str]]:
        """What each cluster the movers leave keeps."""
        return [members_of[cluster].difference(movers) for cluster in {cluster_of(m) for m in movers}]

    def cut_loose(rest: set[str]) -> bool:
        """Whether what stays is one id alone, or in pieces."""
        if len(rest) < 2:
            return bool(rest)
        start = next(iter(rest))
        reached, todo = {start}, [start]
        while todo:
            for other in neighbours[todo.pop()] & (rest - reached):
                reached.add(other)
                todo.append(other)
        return reached != rest

    def move(curies: list[str], target: int) -> None:
        for curie in curies:
            members_of[cluster_of(curie)].discard(curie)
            members_of[target].add(curie)
            moved_to[curie] = target

    structures_restored = 0
    for clique in structure_cliques:
        key = next(m for m in clique if _is_inchikey(m))
        home = cluster_of(key)
        away = [m for m in clique if cluster_of(m) != home]
        if 2 * len(away) <= len(clique):
            continue
        skeleton = _connectivity(key)
        with_its_stereoisomer = all(
            any(_is_inchikey(other) and _connectivity(other) == skeleton for other in members_of[cluster_of(m)])
            for m in away
        )
        if not with_its_stereoisomer or violates([*members_of[home], *away]):
            continue
        # What stays may belong with the movers: ids they tied there, or a SMILES no guardrail weighs against them.
        if any(cut_loose(rest) or any(map(_is_smiles, rest)) for rest in what_stays(away)):
            continue
        move(away, home)
        structures_restored += 1

    def vetted(curie: str, names: set[str], branches: frozenset[str]) -> bool:
        return (
            _prefix(curie) in RETYPABLE_PREFIXES
            and curie not in in_gene_protein_clique
            and not any(map(gene_or_protein_record, members_of[cluster_of(curie)]))
            and not names.isdisjoint(name_keys(label_of(curie), branches))
        )

    for proteins, links in protein_cliques:
        homes = {cluster_of(m) for m in links}
        if len(homes) > 1:
            continue
        home = homes.pop()
        link = min(links)
        branches = info(link).branches
        names = set().union(*(chemical_names[m] for m in links))
        if not all(vetted(m, names, branches) for m in proteins):
            continue
        retype = dict.fromkeys(proteins, branches)
        movers = [m for m in proteins if cluster_of(m) != home]
        # Only a concept alone, or with its clique only, moves: an id beside it is tied to it, not to the chemical.
        if any(what_stays(movers)) or violates([*members_of[home], *movers], retype):
            continue
        move(movers, home)
        retyped_branches.update(retype)
        repairs.retyped.update(dict.fromkeys(proteins, categories[link]))

    for cluster, members in members_of.items():
        if members != set(before[cluster]):
            repairs.replaced.append(before[cluster])
            if members:
                repairs.remade.append(sorted(members))
    logging.info(
        "entity_resolution: chemical Babel cliques: %d moved to their InChIKey, %d protein-typed ids typed as the "
        "chemical they name; %d clusters remade as %d",
        structures_restored,
        len(repairs.retyped),
        len(repairs.replaced),
        len(repairs.remade),
    )
    return repairs
