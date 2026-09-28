"""Memory-bounded entity resolution over the full harmonized graph.

Produces the two artifacts the build needs:

* the integrated **canonical nodes** file (one node per cluster), and
* a ``node_id -> representative_curie`` map for resolving edge endpoints,

replacing ``integrate.py``'s legacy equivalency-trusting node merge.

The build must stay under ~48 GB RAM at ~30M CURIEs, so the two structures that
would otherwise dominate memory are kept **out of core**:

* **evidence accumulation** — all evidence is streamed to a temp file and
  combined with an external ``sort`` (the pattern ``integrate_edges`` already
  uses), so peak memory is one CURIE-pair's evidence, not the whole graph;
* **node materialization** — harmonized nodes are streamed, tagged with their
  cluster, external-sorted by cluster, and reconciled one cluster at a time.

The numeric core (factorization, connected components) uses numpy/scipy, which
are both far faster and far more memory-compact than Python dicts/objects:
CURIE strings become int codes, components come from ``scipy`` csgraph, and label
propagation runs per non-trivial component so peak clustering memory is bounded by
the largest component rather than the whole graph.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import jsonlines
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from kraken.entity_resolution.babel_outliers import (
    compatible_name_pairs,
    drop_corroborated,
    families_signature,
    find_babel_outliers,
    log_babel_outliers,
)
from kraken.entity_resolution.clustering import DEFAULT_SEED, label_propagation
from kraken.entity_resolution.debug_db import DebugDbWriter, debug_db_path
from kraken.entity_resolution.families import ALL_FAMILIES, GENE_PROTEIN_FAMILY, BranchFamilies
from kraken.entity_resolution.gene_protein_cohesion import rejoin_split_gene_protein_cliques
from kraken.entity_resolution.guardrails import (
    GuardrailConfig,
    NodeInfo,
    cluster_violations,
    enforce_cluster,
    ids_per_cluster_histogram,
    log_one_id_repairs,
    log_oversized_clusters,
)
from kraken.entity_resolution.id_facts import IdFacts, IdFactsStore
from kraken.entity_resolution.match_graph import alias_evidence, clique_evidence, match_predicate_evidence
from kraken.entity_resolution.materialize import PrefixRanking, materialize_cluster
from kraken.entity_resolution.name_sim import DEFAULT_STOPLIST, is_droppable, name_keys, normalize_name
from kraken.entity_resolution.prefix_backups import infer_category, infer_taxon
from kraken.entity_resolution.uncanonicalize import (
    ALIAS_EVIDENCE_SOURCES,
    CANONICALIZED_AGGREGATOR_SOURCES,
    is_coarser_than_canonical,
    original_alias_pairs,
)
from kraken.entity_resolution.uncanonicalize import (
    original_endpoints as _original_endpoints,
)
from kraken.entity_resolution.weights import NAME_SIMILARITY_GROUP, ERWeights
from kraken.utils.constants import (
    BABEL_RELATION_ATTRIBUTE,
    EDGE_ATTRIBUTES,
    EDGE_OBJECT,
    EDGE_PREDICATE,
    EDGE_PRIMARY_KS,
    EDGE_SUBJECT,
    EXACT_MATCH_PREDICATES,
    GENE_PROTEIN_CONFLATION_RELATION,
    NODE_CATEGORIES,
    NODE_EQUIVALENT_IDS,
    NODE_ID,
    NODE_NAME,
    NODE_PROVIDED_BY,
    NODE_TAXON,
    SAME_AS_PREDICATE,
)
from kraken.utils.kg_io import remove_file, stream_edges_from_jsonl, stream_nodes_from_jsonl

SEP = "\t"

# Set in an id's ``seeds`` mask when it was seeded ONLY by an alias (not by any node or equivalency list), so the
# edge pass can tell "already in the graph" from "introduced by an earlier alias". Far above any source's bit.
ALIAS_SEED_BIT = 1 << 62


def _evidence_row(a: str, b: str, group: str, weight: float, *, kind: str) -> str:
    """One evidence line. ``kind`` (e.g. "equiv:kg2", "alias:robokop", "babel", "name_sim") is not used to weigh
    anything -- ``group`` does that -- but survives to the oversized-cluster diagnostics, which is the only way to
    tell a Babel clique from an aggregator list once both are in the babel_derived group."""
    return f"{a}{SEP}{b}{SEP}{group}{SEP}{weight}{SEP}{kind}\n"


def _external_sort(input_path: Path, output_path: Path, key_args: list[str], temp_dir: Path) -> None:
    """Byte-ordered (deterministic) external sort with spill files on ``temp_dir``."""
    subprocess.run(
        ["sort", "-t", SEP, *key_args, "-T", str(temp_dir), "-o", str(output_path), str(input_path)],
        check=True,
        env={**os.environ, "LC_ALL": "C"},
    )


# --------------------------------------------------------------------------------------
# Stage 1: stream harmonized data -> evidence file, name file, per-CURIE guardrail facts
# --------------------------------------------------------------------------------------


def _aggregator_name_is_matchable(curie: str, branches: frozenset[str], taxon: str | None) -> bool:
    """Whether an aggregator's name for an id Babel doesn't know may be matched on: anything except a gene or
    protein with no taxon, whose name is a symbol that repeats across species (see the call site)."""
    if branches is ALL_FAMILIES or GENE_PROTEIN_FAMILY not in branches:
        return True
    return bool(taxon or infer_taxon(curie))


def _stage1_write_evidence_and_facts(
    config,
    weights: ERWeights,
    families: BranchFamilies,
    evidence_path: Path,
    names_path: Path,
    source_bits: dict[str, int],
    facts: IdFactsStore,
    *,
    deferred_path: Path | None = None,
    cliques_path: Path | None = None,
) -> tuple[dict[str, frozenset[str]], dict[str, str], set[str], dict[str, int], dict[str, set[str]]]:
    """One streaming pass over harmonized nodes+edges. Writes equivalency-clique
    (native sources only) and match-predicate evidence, writes
    ``name_key<TAB>curie<TAB>families`` rows (one per key; see ``name_keys``) for name similarity, and returns:

    * ``inherited_cats``: curie -> candidate categories, propagated from every node
      whose categories are cleanly SINGLE-family onto each id in that node's
      equivalency list. This is how the vast majority of ids (which only ever appear
      as equiv-list members, never as a primary node) get a category. Conflated
      multi-family nodes are deliberately skipped so their mis-typing doesn't spread.
    * per-CURIE harmonized taxon, the set of harmonized node primary ids, and
      ``seeds`` (every curie seen -- node ids + equiv members -- each of which becomes a node).

    Category strings are interned so the facts dict stays compact (there are only
    ~150 distinct Biolink categories).

    Babel's nodes are also recorded into ``facts`` -- each id's own name, category and taxon -- which later stages
    treat as the source of truth for an id's category and taxon.

    For finding Babel clique outliers (see ``babel_outliers``), the aggregator evidence Babel's say overrules goes to
    ``deferred_path`` instead of the evidence file, and every Babel clique member to ``cliques_path``.
    """
    inherited_cats: dict[str, frozenset[str]] = {}
    node_taxon: dict[str, str] = {}
    node_ids: set[str] = set()
    # seeds: every id we've seen -> a bitmask of the sources that provided it (a node
    # id or an equiv-list member) -- the provenance for retaining EVERY id as a node (singleton if it never merges).
    seeds: dict[str, int] = {}

    # Both maps above hold tens of millions of entries but only a few hundred distinct VALUES (category
    # combinations, source combinations), so every value is interned and shared. A fresh set() per id cost
    # ~200 bytes each, and a fresh int per mask another ~30 -- gigabytes at Babel's scale, which is what ran
    # the 39.6M-id build out of memory.
    shared_category_sets: dict[frozenset[str], frozenset[str]] = {}
    shared_masks: dict[int, int] = {}

    def inherit(curie: str, categories: frozenset[str]) -> None:
        current = inherited_cats.get(curie)
        merged = categories if current is None else current | categories
        inherited_cats[curie] = shared_category_sets.setdefault(merged, merged)

    def mark(curie: str, bits: int) -> int:
        mask = seeds.get(curie, 0) | bits
        mask = shared_masks.setdefault(mask, mask)
        seeds[curie] = mask
        return mask

    # per-source infores provided_by (so a retained bare id can carry real provenance).
    source_provided_by: dict[str, set[str]] = defaultdict(set)

    def babel_categories(curie: str) -> frozenset[str] | None:
        """What Babel's own node for ``curie`` contributes to inheritance -- read from ``facts`` on demand, since
        the node pass doesn't store it per id. Same single-family rule as any other source's node."""
        info = facts.get(curie)
        if info is None or not info.categories:
            return None
        branches = families.branches(info.categories)
        if branches is ALL_FAMILIES or len(branches) != 1:
            return None
        categories = frozenset(sys.intern(c) for c in info.categories)
        return shared_category_sets.setdefault(categories, categories)

    def record_aliases(edge: dict, source: str, bit: int, ev) -> None:
        """Attach the ids an aggregator STARTED from that exist nowhere else in the graph.

        An aggregator edge stores "I resolved X to Y" (see ``original_alias_pairs``). When X appears in NO
        source's nodes or equivalency lists, that pairing is the only thing that can place it, so X is seeded
        like an equiv-list member (Babel name/category/taxon, real provenance, a node of its own if it
        never merges) and the X/Y pairing is written as evidence weighted like a two-id list from that
        aggregator -- merge strength, so X joins Y's cluster. This is how robokop's gtex edges, which originate
        at ``HGVS:`` ids no node carries, reach their ``CAID:`` variants.

        When X IS already in the graph, the alias is skipped entirely: some source's list or a Babel
        clique already places X, and an alias could only repeat that or override it. Every bad merge aliases
        caused came from that case -- thousands of ids fanning into one canonical, or originals paired with
        the wrong endpoint -- because only there can an alias join two entities that already exist.

        "Already in the graph" is decided after EVERY source's nodes have been read (the node pass runs to
        completion before any edges), so it doesn't depend on the order sources are processed in.
        """
        for original, canonical in original_alias_pairs(edge, source):
            mask = seeds.get(original, 0)
            if mask and not mask & ALIAS_SEED_BIT:
                continue  # already placed by some source's nodes or lists
            mark(original, bit | ALIAS_SEED_BIT)
            # The original inherits the canonical id's categories the same way an equiv-list member
            # does (the node pass). Without it, an id from a vocabulary Babel doesn't
            # know (HGVS, CAID) would fall through to NamedThing -- a guardrail wildcard -- so the
            # branch guardrail would go inert on exactly the ids this is introducing.
            canonical_cats = inherited_cats.get(canonical)
            babel_cats = babel_categories(canonical)  # Babel's nodes aren't in inherited_cats (see the node pass)
            if babel_cats:
                canonical_cats = babel_cats if canonical_cats is None else canonical_cats | babel_cats
            if canonical_cats:
                inherit(original, canonical_cats)
            if is_coarser_than_canonical(original, canonical):
                # e.g. a bare rsid (position) stored on one CAID (allele): seeded and typed above so the edge
                # can remap onto the position it was asserted about, but never merged into that one allele.
                continue
            ev_alias = alias_evidence(original, canonical, source, weights)
            if ev_alias is not None:
                ev.write(_evidence_row(*ev_alias, kind=f"alias:{source}"))

    # Babel FIRST, so every other source's lists can be tested against what Babel knows (see below); the rest in
    # a fixed order so a build is reproducible.
    sources = sorted(config.all_harmonized_paths_resolved.items(), key=lambda kv: (kv[0] != "babel", kv[0]))
    babel_facts: list[tuple[str, IdFacts]] = []
    dropped_known_pairs = 0
    with (
        open(evidence_path, "w") as ev,
        open(names_path, "w") as nm,
        open(deferred_path or os.devnull, "w") as deferred,
        open(cliques_path or os.devnull, "w") as cliques,
    ):
        # Pass 1: every source's NODES (ids, equivalency lists, categories, taxa, names).
        for source, (nodes_path, _edges_path) in sources:
            bit = 1 << source_bits[source]
            if Path(nodes_path).exists():
                for node in stream_nodes_from_jsonl(Path(nodes_path)):
                    node_id = node.get(NODE_ID)
                    if not node_id:
                        continue
                    if source == "babel":
                        babel_facts.append((node_id, _babel_id_facts(node)))
                        if len(babel_facts) >= BABEL_FACTS_BATCH:
                            facts.record(babel_facts)
                            babel_facts.clear()
                    elif babel_facts:  # flush before the first non-Babel source, so ``facts`` is complete
                        facts.record(babel_facts)
                        babel_facts.clear()
                    node_ids.add(node_id)
                    mark(node_id, bit)
                    source_provided_by[source].update(node.get(NODE_PROVIDED_BY) or ())
                    equiv_ids = node.get(NODE_EQUIVALENT_IDS) or []
                    for equiv_id in equiv_ids:
                        mark(equiv_id, bit)
                    # Equivalency-clique evidence from EVERY source, weighted per source:
                    # native curated lists are strong (>=tau, merge on their own). An aggregator's list is a
                    # STAR from this node and PREFIX-CAPPED (see match_graph.clique_evidence and
                    # ERWeights.max_ids_per_prefix), and all of them share the "babel_derived" source group so
                    # their Babel echo counts once.
                    #
                    # An aggregator pair BABEL KNOWS BOTH SIDES OF contributes nothing: either Babel already
                    # asserts it (6.3M pairs -- redundant), or Babel deliberately keeps the two apart (564k
                    # pairs), and there we trust Babel. The aggregators canonicalize with drug/chemical
                    # conflation ON, so those 564k are overwhelmingly a compound fused with its salts, its
                    # branded products and its combination products -- metformin with metformin hydrochloride
                    # and Jentadueto, glucose with every stereoform. What an aggregator uniquely offers is the
                    # ids Babel has never heard of (9.7M pairs), and those still count. Nothing is lost from the
                    # graph: integration keeps every cross-cluster assertion as a close_match edge.
                    aggregator_list = source in weights.aggregator_list_sources
                    head_known = aggregator_list and facts.knows(node_id)
                    for a, b, group, weight in clique_evidence(equiv_ids, source, weights, head=node_id):
                        row = _evidence_row(a, b, group, weight, kind=f"equiv:{source}")
                        if aggregator_list:
                            other = b if a == node_id else a
                            if head_known and other != node_id and facts.knows(other):
                                dropped_known_pairs += 1
                                deferred.write(row)
                                continue
                        ev.write(row)
                    # Category inheritance. Not for Babel: its nodes carry each id's own type, which later stages
                    # read from ``facts`` ahead of anything inherited, so recording it again here would only
                    # cost memory -- one entry for every Babel id.
                    cats = node.get(NODE_CATEGORIES) or []
                    branches = families.branches(cats) if cats else ALL_FAMILIES
                    if source != "babel" and cats and branches is not ALL_FAMILIES and len(branches) == 1:
                        interned = frozenset(sys.intern(c) for c in cats)
                        inherit(node_id, interned)
                        for equiv_id in equiv_ids:
                            inherit(equiv_id, interned)
                    taxon = node.get(NODE_TAXON)
                    if taxon:
                        node_taxon[node_id] = sys.intern(taxon)
                    # Name rows for name-similarity are PER-ID. A NATIVE source names
                    # its own id (attributable). A CANONICALIZED aggregator node's name
                    # is the clique's PREFERRED label -- not reliably the canonical id's
                    # own name -- so we do NOT emit it; those ids are named per-id by
                    # their own Babel node instead.
                    #
                    # UNLESS BABEL HAS NEVER HEARD OF THE ID, when that reasoning has nothing to stand on: no
                    # Babel node will ever name it, so without this the id carries NO name into the match graph
                    # at all and can never match anything, however exactly its name agrees with another node's.
                    # KEGG:05012 "Parkinson disease pathway" sat alone beside PANTHER.PATHWAY:P00049, same name,
                    # same type, no evidence between them. 11.0M ids were in that state in the 2.3.0 build, and a
                    # sampled ~21,300 of the name matches it cost them were matches the guardrails would have
                    # allowed (GO:0061769 = GO:0034317 "nicotinate riboside kinase activity",
                    # HANCESTRO:0314 = OBO:HANCESTRO_0314, MESH:D002482 = KEGG.COMPOUND:C00760 "Cellulose").
                    # Except a gene or protein with no taxon, and an id Babel doesn't know rarely has one: its
                    # name is a SYMBOL that repeats across species, and the taxon guardrail -- the only thing
                    # that keeps a cow's TNF out of a human's -- goes wildcard without one. That is 27% of those
                    # matches, and in a sample of 172 not one had a taxon on both sides.
                    if source not in CANONICALIZED_AGGREGATOR_SOURCES or (
                        not facts.knows(node_id) and _aggregator_name_is_matchable(node_id, branches, taxon)
                    ):
                        name = node.get(NODE_NAME)
                        words = normalize_name(name)
                        if not is_droppable(words, min_length=weights.min_name_length, stoplist=DEFAULT_STOPLIST):
                            signature = families_signature(branches)
                            for key in name_keys(name, branches):
                                nm.write(f"{key}{SEP}{node_id}{SEP}{signature}\n")
        if babel_facts:
            facts.record(babel_facts)
            babel_facts.clear()
        if dropped_known_pairs:
            logging.info(
                "entity_resolution: ignored %d aggregator equivalency pairs Babel knows both ids of "
                "(Babel decides those; they survive as close_match edges)",
                dropped_known_pairs,
            )

        # Pass 2: every source's EDGES -- aliases (which need pass 1 complete) and match predicates.
        # Match-predicate (close/exact/same_as) edges are match-graph evidence, but
        # only on their ORIGINAL endpoints (see _original_endpoints): a Babel-
        # canonicalized endpoint would just re-import Babel's clustering. KG2 also
        # needs its close_match down-weighted where subclass edges co-occur, so it
        # gets a dedicated two-phase writer.
        dropped_known_matches = Counter()
        for source, (_nodes_path, edges_path) in sources:
            bit = 1 << source_bits[source]
            if not Path(edges_path).exists():
                continue
            if source == "kg2":
                _write_kg2_match_evidence(Path(edges_path), weights, ev, facts, dropped_known_matches, deferred)
                continue
            if source == "babel":
                _write_babel_evidence(Path(edges_path), weights, ev, cliques)
                continue
            takes_aliases = source in ALIAS_EVIDENCE_SOURCES
            for edge in stream_edges_from_jsonl(Path(edges_path)):
                if takes_aliases:
                    record_aliases(edge, source, bit, ev)
                predicate = edge.get(EDGE_PREDICATE, "")
                if weights.predicate_weight(predicate) is None:
                    continue  # not a usable match predicate; skip before un-canon work
                endpoints = _original_endpoints(edge, source)
                if endpoints is None:
                    continue  # canonicalized aggregator without a known un-canonicalizer
                primary_ks = _primary_ks(edge)
                for subject, object_ in endpoints:
                    ev_edge = match_predicate_evidence(
                        subject, object_, predicate, source, weights, primary_ks=primary_ks
                    )
                    if ev_edge is None:
                        continue
                    row = _evidence_row(*ev_edge, kind=f"match:{source}")
                    if _babel_decides(subject, object_, predicate, source, weights, facts):
                        dropped_known_matches[source] += 1
                        deferred.write(row)
                    else:
                        ev.write(row)
        for source, count in sorted(dropped_known_matches.items()):
            logging.info(
                "entity_resolution: ignored %d %s same_as/exact_match pairs Babel knows both ids of "
                "(Babel decides those; they survive as close_match edges)",
                count,
                source,
            )
    return inherited_cats, node_taxon, node_ids, seeds, source_provided_by


def _babel_decides(a: str, b: str, predicate: str, source: str, weights: ERWeights, facts: IdFactsStore) -> bool:
    """True for an aggregator's same_as / exact_match between two ids Babel knows: like the aggregators' equivalence
    lists, it is Babel's call. Of kg2's 41k such pairs, 9.7k join ids Babel keeps in different cliques, and 4.9k of
    those merged in 2.1.1 at full weight: a GO process with one Reactome reaction, a drug's salt or prodrug with the
    drug (salsalate / salicylic acid), one enantiomer with the other ((S)- / (R)-warfarin)."""
    return (
        predicate in EXACT_MATCH_PREDICATES
        and source in weights.aggregator_list_sources
        and facts.knows(a)
        and facts.knows(b)
    )


def _primary_ks(edge: dict) -> str | None:
    """An edge's primary knowledge source, as a single id (harmonized edges may carry a list)."""
    primary_ks = edge.get(EDGE_PRIMARY_KS)
    if isinstance(primary_ks, list):
        return primary_ks[0] if primary_ks else None
    return primary_ks or None


# How many Babel per-id facts to buffer between sqlite writes.
BABEL_FACTS_BATCH = 100_000


def _babel_id_facts(node: dict) -> IdFacts:
    """A Babel node's own facts."""
    taxon = node.get(NODE_TAXON)
    return IdFacts(
        label=node.get(NODE_NAME),
        categories=tuple(node.get(NODE_CATEGORIES) or ()),
        taxa=(taxon,) if taxon else (),
    )


def _babel_relation(edge: dict) -> str | None:
    """The Babel relation recorded on one of its edges (see BABEL_RELATION_ATTRIBUTE), if any."""
    for attributes in (edge.get(EDGE_ATTRIBUTES) or {}).values():
        if isinstance(attributes, dict) and BABEL_RELATION_ATTRIBUTE in attributes:
            return attributes[BABEL_RELATION_ATTRIBUTE]
    return None


def _conflated_hubs(edges_path: Path) -> dict[str, str]:
    """``clique hub -> the hub the conflated group is anchored on``, from Babel's gene/protein conflation edges.

    Babel writes a gene's clique and its protein's clique separately, then one edge between their preferred ids.
    Treated literally that is two cliques joined by a single bridge -- and a bridge is exactly what a clustering
    step cuts: BRCA1 came out as HGNC + ENSG + OMIM in one node and NCBIGene + its 22 UniProtKB proteins in
    another. We conflate gene and protein, so the two cliques are ONE clique; this map is how the two groups find
    each other (union-find over the conflation edges, so a gene conflated with several proteins is one group).

    Only hubs that take part in a conflation are held, so this is a few million entries, not one per Babel id.
    """
    parent: dict[str, str] = {}

    def find(hub: str) -> str:
        parent.setdefault(hub, hub)
        while parent[hub] != hub:
            parent[hub] = parent[parent[hub]]
            hub = parent[hub]
        return hub

    for edge in stream_edges_from_jsonl(edges_path):
        if _babel_relation(edge) != GENE_PROTEIN_CONFLATION_RELATION:
            continue
        gene, protein = find(edge.get(EDGE_SUBJECT)), find(edge.get(EDGE_OBJECT))
        if gene != protein:
            parent[protein] = gene  # anchor the conflated group on the gene's hub
    return {hub: find(hub) for hub in parent}


def _write_babel_evidence(edges_path: Path, weights: ERWeights, ev, cliques=None) -> None:
    """Equivalence evidence from Babel's same_as edges: each clique, and each gene/protein conflation MERGED INTO
    the gene's clique (see ``_conflated_hubs``), weighted as ONE list with the hub as head -- so a group within
    ``clique_cap`` is a full clique. Babel's other edges -- its drug/chemical relations -- are deliberately not
    evidence: that conflation is off. Each clique's members also go to ``cliques``, if given, as
    ``id<TAB>hub<TAB>clique size`` rows.

    The harmonizer writes each hub's edges consecutively, so an unconflated clique is one run of edges sharing a
    subject and is emitted as it streams. A conflated group's two cliques are written in different places (one per
    compendium), so those are accumulated by anchor hub and emitted at the end; memory is bounded by the
    gene/protein ids that take part in a conflation, not by Babel.
    """
    anchor_of = _conflated_hubs(edges_path)
    logging.info("entity_resolution: %d Babel clique hubs are gene/protein conflations", len(anchor_of))
    conflated: dict[str, set[str]] = defaultdict(set)
    current: str | None = None
    members: list[str] = []

    def record_clique(hub: str, clique: list[str]) -> None:
        if cliques is not None:
            unique = set(clique)
            cliques.writelines(f"{member}{SEP}{hub}{SEP}{len(unique)}\n" for member in unique)

    def flush() -> None:
        if current is None or not members:
            return
        anchor = anchor_of.get(current)
        if anchor is not None:  # part of a conflated group: hold it until every clique of the group is in
            conflated[anchor].update(members)
            conflated[anchor].add(current)
            return
        record_clique(current, [current, *members])
        for evidence in clique_evidence([current, *members], "babel", weights, head=current):
            ev.write(_evidence_row(*evidence, kind="babel"))

    for edge in stream_edges_from_jsonl(edges_path):
        if edge.get(EDGE_PREDICATE) != SAME_AS_PREDICATE:
            continue
        subject = edge.get(EDGE_SUBJECT)
        if subject != current:
            flush()
            current, members = subject, []
        members.append(edge.get(EDGE_OBJECT))
    flush()

    for anchor, group in conflated.items():
        record_clique(anchor, [anchor, *group])
        for evidence in clique_evidence([anchor, *sorted(group)], "babel", weights, head=anchor):
            ev.write(_evidence_row(*evidence, kind="babel:gene_protein"))


# subclass_of / superclass_of between the same pair signals the co-occurring close_match
# is a mislabeled hierarchical relation, not equivalence, so we down-weight it.
SUBCLASS_PREDICATES: frozenset[str] = frozenset({"biolink:subclass_of", "biolink:superclass_of"})


def _subclass_penalized_weight(base_weight: float, hierarchical_count: int, decay: float) -> float:
    """Down-weight a close_match by ``decay`` for each co-occurring hierarchical edge
    (more hierarchical evidence -> weaker close_match)."""
    return base_weight * (decay**hierarchical_count)


def _write_kg2_match_evidence(
    edges_path: Path, weights: ERWeights, ev, facts: IdFactsStore, dropped_known_matches: Counter, deferred=None
) -> None:
    """Emit KG2 match-predicate evidence on un-canonicalized endpoints, down-weighting
    each close_match by how many subclass/superclass edges the same original pair has,
    and setting aside the same_as / exact_match pairs Babel decides (see ``_babel_decides``) -- to ``deferred``, if
    given, in case one turns out to be a Babel clique outlier.

    Two-phase over KG2's edges (option (a)): first count hierarchical edges per original
    pair and buffer the match pairs (bounded by KG2's edge count, not the whole graph),
    then emit each with its penalty applied.

    KG2's originals contribute no alias evidence (see uncanonicalize.ALIAS_EVIDENCE_SOURCES).
    """
    hierarchical_counts: dict[tuple[str, str], int] = defaultdict(int)
    # (a, b, predicate, primary_ks, whether Babel decides it), a <= b
    match_pairs: list[tuple[str, str, str, str | None, bool]] = []
    for edge in stream_edges_from_jsonl(edges_path):
        predicate = edge.get(EDGE_PREDICATE, "")
        is_hierarchical = predicate in SUBCLASS_PREDICATES
        if not is_hierarchical and weights.predicate_weight(predicate) is None:
            continue
        for subject, object_ in _original_endpoints(edge, "kg2") or []:
            if subject == object_:
                continue
            a, b = (subject, object_) if subject <= object_ else (object_, subject)
            if is_hierarchical:
                hierarchical_counts[(a, b)] += 1
            else:
                babel_decides = _babel_decides(a, b, predicate, "kg2", weights, facts)
                if babel_decides:
                    dropped_known_matches["kg2"] += 1
                match_pairs.append((a, b, predicate, _primary_ks(edge), babel_decides))

    for a, b, predicate, primary_ks, babel_decides in match_pairs:
        base = weights.predicate_weight(predicate)
        weight = _subclass_penalized_weight(base, hierarchical_counts.get((a, b), 0), weights.subclass_penalty_decay)
        # Per-primary-KS group, so parallel close_matches on this pair from different KSes sum.
        group = weights.predicate_group("kg2", primary_ks)
        row = _evidence_row(a, b, group, weight, kind="match:kg2")
        if not babel_decides:
            ev.write(row)
        elif deferred is not None:
            deferred.write(row)


def _stage1b_append_name_similarity(
    names_path: Path,
    evidence_path: Path,
    weights: ERWeights,
    temp_dir: Path,
    *,
    guardrail_config: GuardrailConfig | None = None,
    name_pairs_path: Path | None = None,
) -> None:
    """Group CURIEs by normalized name (external sort) and append name-similarity
    clique evidence for each group within the size cap. Bounded memory: one name
    group at a time. The pairs the pairwise guardrails allow also go to ``name_pairs_path``, if given, as
    ``a<TAB>b<TAB>families(a)<TAB>families(b)`` rows, for finding Babel clique outliers."""
    sorted_names = temp_dir / "er_s1_names_sorted.tmp"
    try:
        _external_sort(names_path, sorted_names, ["-k1,1"], temp_dir)
        with (
            open(sorted_names) as fin,
            open(evidence_path, "a") as ev,
            open(name_pairs_path or os.devnull, "w") as name_pairs,
        ):
            current = None
            group: dict[str, str] = {}  # curie -> its families signature

            def flush() -> None:
                unique_ids = sorted(group)
                if not (2 <= len(unique_ids) <= weights.name_group_cap):
                    return
                w = weights.name_similarity_weight
                for i in range(len(unique_ids)):
                    for j in range(i + 1, len(unique_ids)):
                        a, b = unique_ids[i], unique_ids[j]
                        ev.write(_evidence_row(a, b, NAME_SIMILARITY_GROUP, w, kind="name_sim"))
                if name_pairs_path is not None and guardrail_config is not None:
                    for a, b in compatible_name_pairs(group, guardrail_config):
                        name_pairs.write(f"{a}{SEP}{b}{SEP}{group[a]}{SEP}{group[b]}\n")

            for line in fin:
                name, curie, signature = line.rstrip("\n").split(SEP)
                if name != current and group:
                    flush()
                    group = {}
                current = name
                group.setdefault(curie, signature)
            if group:
                flush()
    finally:
        remove_file(sorted_names)


def _restore_deferred_evidence(deferred_path: Path, evidence_path: Path, ids: set[str] | dict) -> int:
    """Append the set-aside aggregator evidence (see ``deferred_path`` in stage 1) that touches any of ``ids``."""
    restored = 0
    with open(deferred_path) as fin, open(evidence_path, "a") as ev:
        for line in fin:
            a, b, _rest = line.split(SEP, 2)
            if a in ids or b in ids:
                ev.write(line)
                restored += 1
    return restored


# --------------------------------------------------------------------------------------
# Stage 2: accumulate evidence -> tau-filtered weighted pairs (external sort)
# --------------------------------------------------------------------------------------


def _stage2_accumulate_pairs(
    evidence_path: Path,
    pairs_path: Path,
    weights: ERWeights,
    temp_dir: Path,
    babel_outliers: set[str] | dict | None = None,
    debug_db: DebugDbWriter | None = None,
) -> int:
    """Combine evidence per CURIE pair (max within source group, sum across)
    and keep pairs meeting tau, leaving out Babel's evidence about its clique outliers (see ``babel_outliers``).
    Every pair with evidence other than a shared Babel clique also goes to ``debug_db``, whether or not it reached
    tau. Returns the number of pairs written."""
    babel_outliers = babel_outliers or set()
    sorted_ev = temp_dir / "er_s2_evidence_sorted.tmp"
    n_pairs = 0
    try:
        _external_sort(evidence_path, sorted_ev, ["-k1,1", "-k2,2", "-k3,3"], temp_dir)
        with open(sorted_ev) as fin, open(pairs_path, "w") as out:
            cur_a = cur_b = None
            group_max: dict[str, float] = {}

            kind_max: dict[str, float] = {}

            def flush() -> int:
                if cur_a is None:
                    return 0
                total = sum(group_max.values())
                if debug_db is not None and not all(kind.startswith("babel") for kind in kind_max):
                    debug_db.add_evidence(cur_a, cur_b, total, total >= weights.tau, kind_max)
                if total >= weights.tau:
                    out.write(f"{cur_a}{SEP}{cur_b}{SEP}{total}\n")
                    return 1
                return 0

            for line in fin:
                a, b, group, weight_s, kind = line.rstrip("\n").split(SEP)
                if babel_outliers and kind.startswith("babel") and (a in babel_outliers or b in babel_outliers):
                    continue
                weight = float(weight_s)
                if a != cur_a or b != cur_b:
                    n_pairs += flush()
                    cur_a, cur_b = a, b
                    group_max = {}
                    kind_max = {}
                prev = group_max.get(group)
                if prev is None or weight > prev:
                    group_max[group] = weight
                if weight > kind_max.get(kind, -1.0):
                    kind_max[kind] = weight
            n_pairs += flush()
    finally:
        remove_file(sorted_ev)
    return n_pairs


# --------------------------------------------------------------------------------------
# Stage 3: components + Leiden + guardrails -> curie -> cluster id
# --------------------------------------------------------------------------------------


def _stage3_cluster(
    pairs_path: Path,
    weights: ERWeights,
    families: BranchFamilies,
    guardrail_config: GuardrailConfig,
    inherited_cats: dict[str, frozenset[str]],
    node_taxon: dict[str, str],
    node_ids: set[str],
    facts: IdFactsStore,
    seed: int,
    cliques_path: Path | None = None,
) -> tuple[dict[str, int], dict[str, dict[int, int]], dict[str, tuple[str, ...]]]:
    """Cluster the weighted pair graph. Returns ``curie -> cluster_id`` for every
    CURIE appearing in a pair, the ids-per-prefix histogram, and ``curie -> single
    intrinsic category`` for every match-graph node.

    Uses int codes (numpy) and scipy connected components; Leiden + guardrails run
    per non-trivial component so peak memory is bounded by the largest component.
    """
    if os.path.getsize(pairs_path) == 0:
        return {}, {}, {}

    df = pd.read_csv(pairs_path, sep=SEP, names=["a", "b", "w"], dtype={"a": str, "b": str, "w": "float32"})
    codes, uniques = pd.factorize(pd.concat([df["a"], df["b"]], ignore_index=True), sort=False)
    n = len(df)
    a_codes = codes[:n].astype(np.int64)
    b_codes = codes[n:].astype(np.int64)
    w = df["w"].to_numpy()
    num_nodes = len(uniques)
    del df

    # Connected components on the undirected graph.
    graph = coo_matrix((w, (a_codes, b_codes)), shape=(num_nodes, num_nodes))
    n_components, labels = connected_components(graph, directed=False)
    logging.info("entity_resolution: %d nodes, %d pairs, %d connected components", num_nodes, n, n_components)

    # Edges grouped by component (both endpoints share a label): sort edge indices
    # by component so each component's edges are a contiguous slice.
    edge_labels = labels[a_codes]
    edge_order = np.argsort(edge_labels, kind="stable")
    edge_labels_sorted = edge_labels[edge_order]
    comp_edge_starts = np.searchsorted(edge_labels_sorted, np.arange(n_components), side="left")
    comp_edge_ends = np.searchsorted(edge_labels_sorted, np.arange(n_components), side="right")

    # Node codes grouped by component.
    node_order = np.argsort(labels, kind="stable")
    labels_sorted = labels[node_order]
    comp_node_starts = np.searchsorted(labels_sorted, np.arange(n_components), side="left")
    comp_node_ends = np.searchsorted(labels_sorted, np.arange(n_components), side="right")

    # One category SET per match-graph node (may hold several biolink categories —
    # we keep all of them, and the node's family is their union), from a chain that
    # never leaves it untyped:
    #   1. Babel's category for this id (source of truth; types each id individually,
    #      so an id a source lists on a conflated node still gets its own type);
    #   2. else the categories inherited from the single-family nodes that list this id
    #      in their equivalency list — the full inherited set, but ONLY if it resolves
    #      to a single family; cross-family disagreement (different aggregators typing
    #      it differently) is ambiguous and falls through;
    #   3. else the prefix->category backup (infer_category) — a LAST-resort guess that
    #      must never pre-empt a real source category, so it runs AFTER inheritance;
    #   4. else NamedThing (a guardrail wildcard).
    # This keeps the branch guardrail from going inert on the huge fraction of ids
    # that only ever appear as equiv-list members, without re-importing Babel's
    # mis-typing (multi-family nodes never propagate in stage 1).
    # One Babel lookup per id as we go, rather than one dict of every match-graph id's facts, and the category
    # tuples interned like stage 1's sets: tens of millions of entries, a few hundred distinct values.
    logging.info("entity_resolution: looking up %d match-graph node categories/taxa", num_nodes)
    mg_categories: dict[str, tuple[str, ...]] = {}
    mg_taxon: dict[str, str] = {}
    shared_category_tuples: dict[tuple[str, ...], tuple[str, ...]] = {}
    for code in range(num_nodes):
        curie = uniques[code]
        norm = facts.get(curie)
        cats = tuple(norm.categories) if norm and norm.categories else ()  # 1. Babel (source of truth)
        if not cats:
            inherited = inherited_cats.get(curie)  # 2. inherited from single-family source equiv-lists
            if inherited:
                branches = families.branches(tuple(inherited))
                if branches is not ALL_FAMILIES and len(branches) == 1:
                    cats = tuple(sorted(inherited))  # single-family inherited category
        if not cats:
            inferred = infer_category(curie)  # 3. prefix backup -- LAST resort, never before source inheritance
            if inferred:
                cats = (inferred,)
        if not cats:
            cats = ("biolink:NamedThing",)  # 4. untyped/ambiguous -> NamedThing (a guardrail wildcard)
        mg_categories[curie] = shared_category_tuples.setdefault(cats, cats)
        # Taxon precedence (mirrors category): 1. Babel (source of truth);
        # 2. the harmonized node's (source) taxon; 3. the single-species prefix backup
        # -- LAST, never before source; else untaxoned (a guardrail wildcard).
        taxa = norm.taxa if norm else ()
        if taxa:
            mg_taxon[curie] = sys.intern(taxa[0])
        elif curie in node_taxon:
            mg_taxon[curie] = node_taxon[curie]
        else:
            inferred_taxon = infer_taxon(curie)
            if inferred_taxon:
                mg_taxon[curie] = inferred_taxon

    def info_provider(curie: str) -> NodeInfo:
        # Stamp the node with its resolved branch-FAMILY set (category -> family done
        # once here), so the guardrails cluster on families directly.
        return NodeInfo(
            curie=curie,
            branches=families.branches(mg_categories.get(curie, ())),
            taxon=mg_taxon.get(curie),
        )

    curie_to_cluster: dict[str, int] = {}
    # The ids-per-prefix histogram and the oversized clusters, tallied as clusters are made rather than by keeping
    # every cluster (and then a sorted copy of every cluster) until the end.
    histogram: dict[str, Counter] = defaultdict(Counter)
    oversized: list[list[str]] = []
    next_cluster_id = 0
    repairs: Counter = Counter()  # one-id repairs, reported once at the end rather than per cluster
    # A cluster that breaks a guardrail is repaired AFTER the component pass, because the repair wants to know
    # which Babel clique each member is in (see greedy_valid_partition) and an id -> hub map for the whole build
    # would not fit in memory. Held with the edges between its members, which is all the repair reads.
    to_repair: list[tuple[list[str], list[tuple[str, str, float]]]] = []

    def record(cluster: list[str]) -> None:
        nonlocal next_cluster_id
        for curie in cluster:
            curie_to_cluster[curie] = next_cluster_id
        for prefix, count in ids_per_cluster_histogram([cluster]).items():
            for ids_of_prefix, clusters in count.items():
                histogram[prefix][ids_of_prefix] += clusters
        if len(cluster) >= guardrail_config.oversized_cluster_log_threshold:
            oversized.append(cluster)
        next_cluster_id += 1

    for comp in range(n_components):
        node_idx = node_order[comp_node_starts[comp] : comp_node_ends[comp]]
        member_curies = [uniques[c] for c in node_idx]
        if len(member_curies) == 1:
            raw_clusters = [member_curies]
        else:
            edge_idx = edge_order[comp_edge_starts[comp] : comp_edge_ends[comp]]
            comp_edges = [(uniques[a_codes[i]], uniques[b_codes[i]], float(w[i])) for i in edge_idx]
            info = {c: info_provider(c) for c in member_curies}
            # Guardrail-aware FORMATION: drop edges between conflicting nodes before
            # clustering, so incompatible things never merge in the first place (this
            # can also split the component for free). cluster_violations on the pair
            # covers all enforced guardrails: branch, one_id, taxon.
            comp_edges = [e for e in comp_edges if not cluster_violations([e[0], e[1]], info, guardrail_config)]
            # Label propagation on the pruned component (no resolution parameter —
            # LP merges what's connected; the guardrails do the splitting).
            raw_clusters = label_propagation(member_curies, comp_edges, seed=seed)
            # Guardrails as the backstop: a cluster that still violates one (a transitive conflict the pairwise
            # prune above cannot see) is held back for the repair pass after this loop, with the edges between
            # its members -- LP has no resolution to raise, so greedy_valid_partition is the splitter.
            checked: list[list[str]] = []
            for cluster in raw_clusters:
                if cluster_violations(sorted(cluster), info, guardrail_config):
                    held = set(cluster)
                    to_repair.append(
                        (cluster, [e for e in comp_edges if e[0] in held and e[1] in held]),
                    )
                else:
                    checked.append(cluster)
            raw_clusters = checked

        for cluster in raw_clusters:
            record(cluster)

    # The guardrail repairs, now that the Babel cliques the violating clusters are made of can be looked up.
    if to_repair:
        wanted = {curie for cluster, _edges in to_repair for curie in cluster}
        clique_of = _clique_membership(cliques_path, wanted, weights.clique_cap) if cliques_path else {}
        logging.info(
            "entity_resolution: repairing %d guardrail-violating clusters (%d ids, %d of them in a Babel clique "
            "the repair keeps whole)",
            len(to_repair),
            len(wanted),
            len(clique_of),
        )
        for cluster, edges in to_repair:
            info = {c: info_provider(c) for c in cluster}
            for part in enforce_cluster(
                sorted(cluster),
                info,
                guardrail_config,
                adjacency=_adjacency(edges),
                splitter=None,
                repairs=repairs,
                clique_of=clique_of,
            ):
                record(part)

    # A Babel gene/protein clique that label propagation tore in two is put back (see gene_protein_cohesion),
    # and what it counted per cluster retallied for the clusters that changed.
    if cliques_path is not None:
        for group in rejoin_split_gene_protein_cliques(curie_to_cluster, cliques_path, info_provider, guardrail_config):
            for cluster in group:
                for prefix, count in ids_per_cluster_histogram([cluster]).items():
                    for ids_of_prefix, clusters in count.items():
                        histogram[prefix][ids_of_prefix] -= clusters
            rejoined = [curie for cluster in group for curie in cluster]
            for prefix, count in ids_per_cluster_histogram([rejoined]).items():
                for ids_of_prefix, clusters in count.items():
                    histogram[prefix][ids_of_prefix] += clusters
            if len(rejoined) >= guardrail_config.oversized_cluster_log_threshold:
                oversized.append(rejoined)

    log_one_id_repairs(repairs, guardrail_config)
    log_oversized_clusters(oversized, guardrail_config)
    return curie_to_cluster, {prefix: dict(counts) for prefix, counts in histogram.items()}, mg_categories


def _clique_membership(cliques_path: Path, wanted: set[str], clique_cap: int) -> dict[str, str]:
    """``id -> Babel clique hub`` for the wanted ids, from the ``id, hub, size`` rows stage 1 wrote.

    Only cliques within ``clique_cap`` are returned: a larger one is emitted as a star from its hub (see
    ``match_graph.clique_evidence``), so its members are connected only through the hub and treating them as one
    group would assert a link the evidence does not have.
    """
    membership: dict[str, str] = {}
    with open(cliques_path) as fin:
        for line in fin:
            curie, hub, size = line.rstrip("\n").split(SEP)
            if curie in wanted and int(size) <= clique_cap:
                membership[curie] = sys.intern(hub)
    return membership


def _adjacency(edges: list[tuple[str, str, float]]) -> dict[str, dict[str, float]]:
    adj: dict[str, dict[str, float]] = defaultdict(dict)
    for a, b, weight in edges:
        adj[a][b] = weight
        adj[b][a] = weight
    return adj


# --------------------------------------------------------------------------------------
# Stage 4: materialize canonical nodes (external group-by) -> node_id -> representative
# --------------------------------------------------------------------------------------


def _stage4_materialize(
    config,
    curie_to_cluster: dict[str, int],
    mg_categories: dict[str, tuple[str, ...]],
    node_ids: set[str],
    seeds: dict[str, int],
    source_provided_by: dict[str, set[str]],
    source_bits: dict[str, int],
    inherited_cats: dict[str, frozenset[str]],
    facts: IdFactsStore,
    ranking: PrefixRanking,
    families: BranchFamilies,
    biolink,
    temp_dir: Path,
    debug_db: DebugDbWriter | None = None,
) -> dict[str, str]:
    """Stream harmonized nodes, group by cluster on disk, reconcile one cluster at a
    time, and write the canonical nodes file (and each id's row of ``debug_db``, if given). Returns
    ``node_id -> representative``.

    EVERY id a source provided is materialized -- merged into its cluster, or emitted
    as its own SINGLETON if it never merged. Bare ids (equiv-list members with no
    harmonized node) become synthetic member dicts carrying their Babel
    name/category/taxon and real provenance (decoded from the ``seeds`` bitmask), so
    nothing a source provided is ever dropped. Peak memory is one cluster's members.
    """
    keyed = temp_dir / "er_s4_nodes_keyed.tmp"
    keyed_sorted = temp_dir / "er_s4_nodes_keyed_sorted.tmp"
    node_id_to_rep: dict[str, str] = {}
    bit_to_source = {bit: src for src, bit in source_bits.items()}

    def provenance(mask: int) -> list[str]:
        provided: set[str] = set()
        for bit, src in bit_to_source.items():
            if mask & (1 << bit):
                provided.update(source_provided_by.get(src, ()) or {src})
        return sorted(provided)

    def label_of(curie: str) -> str | None:
        info = facts.get(curie)
        return info.label if info else None

    def taxon_of(curie: str) -> str | None:
        info = facts.get(curie)
        return info.taxa[0] if info and info.taxa else None

    def node_category(curie: str, source_categories: list[str] | None) -> list[str]:
        """The categories to write on a HARMONIZED node: its intrinsic one, else Babel's, else what its source
        said, else the bare-id chain.

        Only ids that reached the match graph have an intrinsic category, and an id whose evidence never reached
        tau doesn't -- so taking ``mg_categories`` alone silently dropped the type of every unmerged node, half
        the graph: 2.1M NCBITaxon ids Babel calls OrganismTaxon, 1.1M UniProtKB proteins, 169k RefMet small
        molecules, all written out as NamedThing. Babel comes before the source's own list for the same reason
        ``mg_categories`` does: an aggregator's list is conflated, Babel's per-id type is not.
        """
        cats = mg_categories.get(curie)
        if cats:
            return sorted(cats)
        info = facts.get(curie)
        if info and info.categories:
            return sorted(info.categories)
        if source_categories:
            return sorted(source_categories)
        return bare_category(curie)

    def bare_category(curie: str) -> list[str]:
        cats = mg_categories.get(curie)  # computed in stage 3 for match-graph ids
        if cats:
            return sorted(cats)
        info = facts.get(curie)  # isolated id: Babel -> inherited -> prefix backup -> NamedThing
        if info and info.categories:
            return sorted(info.categories)
        inherited = inherited_cats.get(curie)
        if inherited:
            branches = families.branches(tuple(inherited))
            if branches is not ALL_FAMILIES and len(branches) == 1:
                return sorted(inherited)
        inferred = infer_category(curie)  # prefix backup only after source inheritance (see mg_categories)
        if inferred:
            return [inferred]
        return ["biolink:NamedThing"]

    try:
        with jsonlines.open(keyed, "w") as writer:
            for source, (nodes_path, _edges) in sorted(config.all_harmonized_paths_resolved.items()):
                if not Path(nodes_path).exists():
                    continue
                for node in stream_nodes_from_jsonl(Path(nodes_path)):
                    node_id = node.get(NODE_ID)
                    if not node_id:
                        continue
                    if debug_db is not None:
                        node[SOURCE_MARKER] = source  # which source said what, for the debug db; removed at flush
                    # Replace the (possibly conflated) source category list with this id's single intrinsic
                    # category, so the merged node's categories are the union of its members' true types, not
                    # conflation leftovers.
                    node[NODE_CATEGORIES] = node_category(node_id, node.get(NODE_CATEGORIES))
                    # Taxon: Babel wins; else keep the source taxon already on the dict; else
                    # the single-species prefix backup (last resort, never before source).
                    babel_taxon = taxon_of(node_id)
                    if babel_taxon:
                        node[NODE_TAXON] = babel_taxon
                    elif not node.get(NODE_TAXON):
                        inferred_taxon = infer_taxon(node_id)
                        if inferred_taxon:
                            node[NODE_TAXON] = inferred_taxon
                    cid = curie_to_cluster.get(node_id)
                    key = f"c{cid}" if cid is not None else f"s:{node_id}"
                    writer.write([key, node])
            # Synthetic member dicts for EVERY bare id (equiv-list member with no
            # harmonized node): merged into its cluster if it merged, else its own
            # SINGLETON -- so every id a source provided survives to the graph.
            for curie, mask in seeds.items():
                if curie in node_ids:
                    continue  # has a harmonized node (already written above)
                cid = curie_to_cluster.get(curie)
                key = f"c{cid}" if cid is not None else f"s:{curie}"
                synthetic = {
                    NODE_ID: curie,
                    NODE_CATEGORIES: bare_category(curie),
                    NODE_PROVIDED_BY: provenance(mask),
                }
                info = facts.get(curie)
                if info and info.label:
                    synthetic[NODE_NAME] = info.label
                if info and info.taxa:
                    synthetic[NODE_TAXON] = info.taxa[0]
                else:
                    inferred_taxon = infer_taxon(curie)  # bare id: no source taxon; prefix backup last
                    if inferred_taxon:
                        synthetic[NODE_TAXON] = inferred_taxon
                writer.write([key, synthetic])
            # A canonical node's equivalent_ids must be its exact cluster membership, NOT the union of its members'
            # source lists (which can hold ids the guardrails split into other clusters, breaking disjointness).
            # Membership goes through the external sort alongside the nodes, as one marker per id, rather than as
            # an in-memory inverse of curie_to_cluster -- a second copy of the whole match graph, held right when
            # this stage runs out of memory.
            for curie, cid in curie_to_cluster.items():
                writer.write([f"c{cid}", {CLUSTER_MEMBER_MARKER: curie}])
        _external_sort(keyed, keyed_sorted, ["-k1,1"], temp_dir)

        with (
            jsonlines.open(keyed_sorted) as reader,
            jsonlines.open(config.integrated_nodes_path, "w") as out,
        ):
            current_key: str | None = None
            members: list[dict] = []

            def flush(group_key: str, entries: list[dict]) -> None:
                membership = [entry[CLUSTER_MEMBER_MARKER] for entry in entries if CLUSTER_MEMBER_MARKER in entry]
                group = [entry for entry in entries if CLUSTER_MEMBER_MARKER not in entry]
                if not group:
                    return  # a cluster of ids no source provided a node or list entry for: nothing to write
                sources = [entry.pop(SOURCE_MARKER, None) for entry in group]
                node = materialize_cluster(group, ranking, families, label_of=label_of)
                # No drop rule: every id a source provided is kept (bare ids now carry
                # real provenance, so nothing is provenance-less). A bare-only cluster
                # becomes its own node instead of vanishing.
                if biolink is not None and node.get(NODE_CATEGORIES):
                    leaves = biolink.filter_to_leaf_categories(node[NODE_CATEGORIES])
                    if leaves:
                        node[NODE_CATEGORIES] = sorted(leaves)
                # Every node needs a category; an untyped survivor becomes NamedThing.
                if not node.get(NODE_CATEGORIES):
                    node[NODE_CATEGORIES] = ["biolink:NamedThing"]
                # equivalent_ids = exact cluster membership (guarantees disjointness)
                if group_key.startswith("c"):
                    node[NODE_EQUIVALENT_IDS] = sorted(membership)
                else:  # singleton: the node's own id(s)
                    node[NODE_EQUIVALENT_IDS] = sorted({m[NODE_ID] for m in group})
                out.write(node)
                rep = node[NODE_ID]
                for member in group:
                    node_id_to_rep[member[NODE_ID]] = rep
                if debug_db is not None:
                    _record_debug_ids(debug_db, rep, node[NODE_EQUIVALENT_IDS], group, sources)

            for key, node in reader:
                if key != current_key and members:
                    flush(current_key, members)
                    members = []
                current_key = key
                members.append(node)
            if current_key is not None:
                flush(current_key, members)
    finally:
        remove_file(keyed)
        remove_file(keyed_sorted)
    return node_id_to_rep


def _record_debug_babel_cliques(debug_db: DebugDbWriter, cliques_path: Path) -> None:
    """Every Babel clique member (``id<TAB>hub<TAB>size`` rows written in stage 1) into the debug db."""
    with open(cliques_path) as fin:
        for line in fin:
            curie, hub, size = line.rstrip("\n").split(SEP)
            debug_db.add_babel_clique_member(curie, hub, int(size))


def _record_debug_ids(
    debug_db: DebugDbWriter, rep: str, ids: list[str], group: list[dict], sources: list[str | None]
) -> None:
    """One debug-db row per id of a written node, from the source records ER merged it from. Babel's record is the
    id's own name and type; every other source's is what that source called it. A bare id (an equivalence-list
    member no source has a node for) has one synthetic record, named by Babel if Babel knows it."""
    records: dict[str, list[tuple[str | None, dict]]] = defaultdict(list)
    for source, entry in zip(sources, group, strict=True):
        records[entry[NODE_ID]].append((source, entry))
    for curie in ids:
        babel = next((entry for source, entry in records[curie] if source == "babel"), None)
        entries = [entry for _source, entry in records[curie]]
        debug_db.add_id(
            curie,
            rep,
            babel_name=babel.get(NODE_NAME) if babel else None,
            babel_categories=babel.get(NODE_CATEGORIES) or [] if babel else [],
            taxon=next((e[NODE_TAXON] for e in entries if e.get(NODE_TAXON)), None),
            categories=sorted({c for e in entries for c in e.get(NODE_CATEGORIES) or []}),
            provided_by={p for e in entries for p in e.get(NODE_PROVIDED_BY) or []},
            source_names=[
                (source or "(bare id)", entry[NODE_NAME])
                for source, entry in records[curie]
                if entry.get(NODE_NAME) and source != "babel"
            ],
        )


# --------------------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------------------


def _stage_banner(msg: str) -> None:
    """Log a prominent, greppable banner so the numbered ER stages stand out in the very
    long build log (otherwise they're lost among hundreds of thousands of INFO lines)."""
    logging.info("#" * 100)
    logging.info("###  %s", msg)


# Every cluster at least this large gets a breakdown of the evidence that built it, written to
# ``<integrated debug dir>/oversized_clusters.jsonl``. 2.1.1's largest was 1,455 members, so a healthy build writes a
# handful; a regression (the 56,515-member clusters of the first size-aware build) shows what caused it.
OVERSIZED_EVIDENCE_REPORT_MIN_SIZE = 1000
OVERSIZED_REPORT_FILENAME = "oversized_clusters.jsonl"
# Every Babel clique outlier, with the clique its name matches, in ``<integrated debug dir>`` (see babel_outliers).
BABEL_OUTLIERS_REPORT_FILENAME = "babel_clique_outliers.tsv"
# Key of the membership markers stage 4 sorts alongside the node dicts (not a node property).
CLUSTER_MEMBER_MARKER = "__cluster_member__"
# Key stage 4 tags each source record with, for the debug db (removed before the node is built).
SOURCE_MARKER = "__source__"
# Babel's per-id facts for this ER run (see resolve_entities); a temp file, removed when ER finishes.
BABEL_FACTS_FILENAME = "er_babel_facts.tmp.sqlite"
OVERSIZED_REPORT_SAMPLE_IDS = 25


def _report_oversized_cluster_evidence(
    evidence_path: Path, curie_to_cluster: dict[str, int], report_path: Path
) -> None:
    """Write, for every cluster of OVERSIZED_EVIDENCE_REPORT_MIN_SIZE+ members, how many evidence rows of each kind link
    members INSIDE it, plus its id-prefix mix and a sample of ids -- one JSON line per cluster, largest first.

    One pass over the evidence file, and only when such a cluster exists. The kind column (``equiv:<source>``,
    ``alias:<source>``, ``match:<source>``, ``babel``, ``name_sim``) is what separates a Babel clique from an
    aggregator list; the de-correlation group alone can't.
    """
    sizes = Counter(curie_to_cluster.values())
    ranked = [cid for cid, n in sizes.most_common() if n >= OVERSIZED_EVIDENCE_REPORT_MIN_SIZE]
    if not ranked:
        return
    wanted = set(ranked)
    kinds: dict[int, Counter] = {cid: Counter() for cid in ranked}
    prefixes: dict[int, Counter] = {cid: Counter() for cid in ranked}
    samples: dict[int, list[str]] = defaultdict(list)
    for curie, cid in curie_to_cluster.items():
        if cid in wanted:
            prefixes[cid][curie.split(":", 1)[0]] += 1
            if len(samples[cid]) < OVERSIZED_REPORT_SAMPLE_IDS:
                samples[cid].append(curie)
    with open(evidence_path) as evidence:
        for line in evidence:
            a, b, _group, _weight, kind = line.rstrip("\n").split(SEP)
            cid = curie_to_cluster.get(a)
            if cid in wanted and curie_to_cluster.get(b) == cid:
                kinds[cid][kind] += 1
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with jsonlines.open(report_path, "w") as report:
        for cid in ranked:
            report.write(
                {
                    "members": sizes[cid],
                    "evidence_rows_inside_by_kind": dict(kinds[cid].most_common()),
                    "id_prefixes": dict(prefixes[cid].most_common()),
                    "sample_ids": sorted(samples[cid]),
                }
            )
    logging.warning(
        "entity_resolution: %d cluster(s) have %d+ members (largest %d); evidence breakdown written to %s",
        len(ranked),
        OVERSIZED_EVIDENCE_REPORT_MIN_SIZE,
        sizes[ranked[0]],
        report_path,
    )


def _debug_dir(config) -> Path:
    """Where diagnostics go: the build's integrated debug dir (falls back to the integrated dir for ad-hoc configs)."""
    return Path(getattr(config, "integrated_debug_dir", None) or config.integrated_dir)


def resolve_entities(config, biolink) -> dict[str, str]:
    """Run entity resolution end to end (out of core) and write the canonical nodes
    file. Returns ``node_id -> representative_curie`` for edge resolution."""
    weights = ERWeights.load()
    families = BranchFamilies.load()
    ranking = PrefixRanking.load()
    guardrail_config = GuardrailConfig()
    temp_dir = config.integrated_dir
    temp_dir.mkdir(parents=True, exist_ok=True)

    evidence_path = temp_dir / "er_s1_evidence.tmp"
    names_path = temp_dir / "er_s1_names.tmp"
    deferred_path = temp_dir / "er_s1_deferred_evidence.tmp"
    cliques_path = temp_dir / "er_s1_babel_cliques.tmp"
    name_pairs_path = temp_dir / "er_s1b_name_pairs.tmp"
    pairs_path = temp_dir / "er_s2_pairs.tmp"

    # Deterministic source -> bit index, for encoding per-id provenance in ``seeds``.
    source_bits = {src: i for i, src in enumerate(sorted(config.all_harmonized_paths_resolved))}

    # Per-id names/categories/taxa come from Babel's harmonized nodes, recorded in stage 1.
    facts_path = temp_dir / BABEL_FACTS_FILENAME
    remove_file(facts_path)  # rebuilt from this build's Babel nodes every run
    facts = IdFactsStore(facts_path)
    debug_db = DebugDbWriter(debug_db_path(_debug_dir(config), getattr(config, "kraken_version", None)))
    try:
        t = time.perf_counter()
        _stage_banner("ER STAGE 1 -- streaming harmonized nodes/edges -> match evidence, names, guardrail facts")
        inherited_cats, node_taxon, node_ids, seeds, source_provided_by = _stage1_write_evidence_and_facts(
            config,
            weights,
            families,
            evidence_path,
            names_path,
            source_bits,
            facts,
            deferred_path=deferred_path,
            cliques_path=cliques_path,
        )
        _stage_banner(
            f"ER STAGE 1 DONE ({time.perf_counter() - t:.1f}s) -- {len(node_ids)} harmonized node ids, "
            f"{len(seeds)} total ids (incl. equiv-list members)"
        )

        t = time.perf_counter()
        _stage_banner("ER STAGE 1b -- adding name-similarity match evidence, finding Babel clique outliers")
        _stage1b_append_name_similarity(
            names_path,
            evidence_path,
            weights,
            temp_dir,
            guardrail_config=guardrail_config,
            name_pairs_path=name_pairs_path,
        )
        candidates = find_babel_outliers(name_pairs_path, cliques_path, weights.clique_cap, temp_dir)
        outliers = drop_corroborated(candidates, cliques_path, deferred_path)
        logging.info(
            "entity_resolution: %d Babel clique outlier candidates, %d left alone because the aggregators place "
            "them where Babel does",
            len(candidates),
            len(candidates) - len(outliers),
        )
        log_babel_outliers(outliers, _debug_dir(config) / BABEL_OUTLIERS_REPORT_FILENAME)
        if outliers:
            restored = _restore_deferred_evidence(deferred_path, evidence_path, outliers)
            logging.info("entity_resolution: restored %d aggregator evidence rows about Babel outliers", restored)
        _record_debug_babel_cliques(debug_db, cliques_path)
        _stage_banner(f"ER STAGE 1b DONE ({time.perf_counter() - t:.1f}s)")

        t = time.perf_counter()
        _stage_banner("ER STAGE 2 -- accumulating + tau-filtering weighted match pairs")
        n_pairs = _stage2_accumulate_pairs(evidence_path, pairs_path, weights, temp_dir, outliers, debug_db)
        _stage_banner(f"ER STAGE 2 DONE ({time.perf_counter() - t:.1f}s) -- {n_pairs} pairs above tau")

        t = time.perf_counter()
        _stage_banner("ER STAGE 3 -- clustering (connected components -> label propagation -> guardrails)")
        curie_to_cluster, histogram, mg_categories = _stage3_cluster(
            pairs_path,
            weights,
            families,
            guardrail_config,
            inherited_cats,
            node_taxon,
            node_ids,
            facts,
            DEFAULT_SEED,
            cliques_path,
        )
        _stage_banner(
            f"ER STAGE 3 DONE ({time.perf_counter() - t:.1f}s) -- "
            f"{len(curie_to_cluster)} ids -> {len(set(curie_to_cluster.values()))} clusters"
        )
        _log_histogram(histogram)
        _report_oversized_cluster_evidence(
            evidence_path, curie_to_cluster, _debug_dir(config) / OVERSIZED_REPORT_FILENAME
        )

        del node_taxon  # stage 3's alone; stage 4 is where memory is tightest

        t = time.perf_counter()
        _stage_banner(f"ER STAGE 4 -- materializing canonical nodes -> {config.integrated_nodes_path}")
        node_id_to_rep = _stage4_materialize(
            config,
            curie_to_cluster,
            mg_categories,
            node_ids,
            seeds,
            source_provided_by,
            source_bits,
            inherited_cats,
            facts,
            ranking,
            families,
            biolink,
            temp_dir,
            debug_db,
        )
        _stage_banner(
            f"ER STAGE 4 DONE ({time.perf_counter() - t:.1f}s) -- "
            f"{len(node_id_to_rep)} node ids -> {len(set(node_id_to_rep.values()))} canonical nodes"
        )
    finally:
        facts.close()
        debug_db.close()
        remove_file(facts_path)
        for path in (evidence_path, names_path, deferred_path, cliques_path, name_pairs_path, pairs_path):
            remove_file(path)
    logging.info("entity_resolution: debug database (inspect any id's node with it) at %s", debug_db.path)

    _report_eval(curie_to_cluster)
    logging.info("entity_resolution: %d node ids mapped to representatives", len(node_id_to_rep))
    return node_id_to_rep


def _log_histogram(histogram: dict[str, dict[int, int]]) -> None:
    for prefix in ("HGNC", "NCBIGene", "MONDO", "UniProtKB", "RM", "LM"):
        counts = histogram.get(prefix)
        if counts:
            multi = {k: v for k, v in sorted(counts.items()) if k > 1}
            if multi:
                logging.info("entity_resolution ids-per-cluster[%s]: >1 -> %s", prefix, multi)


def _report_eval(curie_to_cluster: dict[str, int]) -> None:
    from kraken.entity_resolution.eval.scorer import DEFAULT_GROUND_TRUTH_PATH, load_gold, score

    if not DEFAULT_GROUND_TRUTH_PATH.exists() or not curie_to_cluster:
        return
    res = score(load_gold(DEFAULT_GROUND_TRUTH_PATH), curie_to_cluster)
    logging.info(
        "entity_resolution eval: precision=%.4f recall=%.4f f1=%.4f (ml covered=%d, cl covered=%d)",
        res.precision,
        res.recall,
        res.f1,
        res.must_link_covered,
        res.cannot_link_covered,
    )
