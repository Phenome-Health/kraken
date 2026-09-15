"""Integration: run the real out-of-core build on real CURIEs from the seed.

For each conflation case we write harmonized nodes for the true entities plus the
weak aggregator clique that caused the merge, run the build, and confirm it splits
them back apart (scored against the ground truth). Exercises the same code path
the actual build ships (build.resolve_entities), including edge pruning.
"""

import json
from pathlib import Path
from types import SimpleNamespace

import jsonlines
import pytest

from kraken.entity_resolution.build import resolve_entities
from kraken.entity_resolution.eval.scorer import load_gold, score

pytest.importorskip("igraph")

# Real Parkinson-disease clique (from the 2.1.0 merged node): a single-entity
# disease clique whose cross-ontology members are carried ONLY by an aggregator's
# equivalency list -- no independent high-weight source. This is the class of
# merge that regressed to singletons when aggregator equivalency sat below tau.
PARKINSON = [
    "MONDO:0005180", "DOID:14330", "UMLS:C0030567", "MESH:D010300",
    "SNOMEDCT:49049000", "NCIT:C26845", "medgen:10590", "ICD9:332",
    "KEGG.DISEASE:05012", "MEDDRA:10061536",
]

REPO_ROOT = Path(__file__).resolve().parents[2]
CFG = REPO_ROOT / "config" / "entity_resolution"
SEED = CFG / "ground_truth_seed.jsonl"


def _category_for(curie: str, group_label: str) -> str:
    prefix = curie.split(":", 1)[0]
    if "Gene / Protein" in group_label:
        return "biolink:Protein" if prefix in {"UniProtKB", "PR"} else "biolink:Gene"
    if "Disease" in group_label:
        return "biolink:Disease"
    if "Chemical" in group_label:
        return "biolink:ChemicalEntity"
    return "biolink:NamedThing"


def _write_source(tmp_path: Path, source: str, records: list[dict], edges: list[dict] = ()) -> tuple[Path, Path]:
    d = tmp_path / "harmonized" / source
    d.mkdir(parents=True, exist_ok=True)
    nodes_path = d / "nodes.jsonl"
    with jsonlines.open(nodes_path, "w") as w:
        w.write_all(records)
    edges_path = d / "edges.jsonl"
    if edges:
        with jsonlines.open(edges_path, "w") as w:
            w.write_all(edges)
    return nodes_path, edges_path


def _babel_node(curie: str, category: str, name: str | None = None, taxon: str | None = None) -> dict:
    node = {"id": curie, "categories": [category], "provided_by": ["infores:sri-node-normalizer"]}
    if name:
        node["name"] = name
    if taxon:
        node["taxon"] = taxon
    return node


def _same_as_star(members: list[str]) -> list[dict]:
    """A Babel clique as the harmonizer writes it: same_as from the preferred (first) id to each other member."""
    return [
        {"subject": members[0], "predicate": "biolink:same_as", "object": member, "primary_knowledge_source": "x"}
        for member in members[1:]
    ]


def _config(tmp_path: Path, harmonized: dict) -> SimpleNamespace:
    integrated = tmp_path / "integrated"
    return SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )


def _build_config_from_seed(tmp_path: Path):
    """Materialize the conflation cases as harmonized sources on disk."""
    gene_nodes: list[dict] = []
    disease_nodes: list[dict] = []
    conflation_nodes: list[dict] = []
    raw = [json.loads(line) for line in open(SEED) if line.strip()]
    for rec in raw:
        if "clusters" not in rec:
            continue
        reps: list[str] = []
        for label, ids in rec["clusters"].items():
            ids = sorted(ids)
            if not ids:
                continue
            reps.append(ids[0])
            bucket = gene_nodes if "Gene" in label else disease_nodes
            source = "ncbigene" if "Gene" in label else "umls"
            for curie in ids:
                bucket.append(
                    {
                        "id": curie,
                        "categories": [_category_for(curie, label)],
                        "equivalent_ids": ids,
                        "provided_by": [source],
                    }
                )
        if len(reps) > 1:  # the weak aggregator conflation spanning true entities
            conflation_nodes.append({"id": reps[0], "categories": [], "equivalent_ids": reps, "provided_by": ["kg2"]})

    harmonized = {
        "ncbigene": _write_source(tmp_path, "ncbigene", gene_nodes),
        "umls": _write_source(tmp_path, "umls", disease_nodes),
        "kg2": _write_source(tmp_path, "kg2", conflation_nodes),
    }
    integrated = tmp_path / "integrated"
    return SimpleNamespace(
        all_harmonized_paths_resolved=harmonized,
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )


def test_build_splits_seed_conflations(tmp_path):
    config = _build_config_from_seed(tmp_path)
    node_id_to_rep = resolve_entities(config, biolink=None)

    gold = [g for g in load_gold(SEED) if g.case.startswith("conflation:")]
    res = score(gold, node_id_to_rep)  # rep acts as the cluster id
    assert res.fp == 0, f"unexpected cross-entity merges: precision={res.precision}"
    assert res.precision == 1.0
    assert res.recall >= 0.99, f"recall too low: {res.recall}"


def test_ace_rtd_specifically_split(tmp_path):
    config = _build_config_from_seed(tmp_path)
    m = resolve_entities(config, biolink=None)
    assert m["NCBIGene:1636"] == m["HGNC:2707"] == m["UniProtKB:P12821"]
    assert m["MONDO:0017609"] == m["orphanet:3033"]
    assert m["NCBIGene:1636"] != m["MONDO:0017609"]


def _aggregator_only_clique_config(tmp_path: Path, members: list[str]):
    """One kg2 node whose equivalency list IS the clique (single Disease family),
    with no other source. Mirrors a disease clique carried only by an aggregator."""
    node = {
        "id": members[0],
        "categories": ["biolink:Disease"],
        "equivalent_ids": sorted(members),
        "provided_by": ["kg2"],
    }
    integrated = tmp_path / "integrated"
    return SimpleNamespace(
        all_harmonized_paths_resolved={"kg2": _write_source(tmp_path, "kg2", [node])},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )


def test_babel_clique_merges_disease(tmp_path):
    """A cross-ontology disease clique merges into one node when Babel supplies the clique (the Parkinson case)."""
    babel = _write_source(
        tmp_path, "babel", [_babel_node(c, "biolink:Disease") for c in PARKINSON], _same_as_star(PARKINSON)
    )
    m = resolve_entities(_config(tmp_path, {"babel": babel}), biolink=None)
    mapped = [c for c in PARKINSON if m.get(c) == m.get(PARKINSON[0])]
    assert mapped == PARKINSON, (
        f"clique did not fully merge: {len(mapped)}/{len(PARKINSON)} merged; "
        f"missing/split={[c for c in PARKINSON if m.get(c) != m.get(PARKINSON[0])]}"
    )


def test_small_aggregator_equiv_list_merges_on_its_own(tmp_path):
    """An aggregator list with no bulk prefix merges without any Babel clique (the kg2 Parkinson clique
    is 10 ids of different prefixes)."""
    config = _aggregator_only_clique_config(tmp_path, PARKINSON)
    m = resolve_entities(config, biolink=None)  # no Babel source -> no cliques
    reps = {m.get(c) for c in PARKINSON}
    assert None not in reps, "an id was dropped"
    assert len(reps) == 1, f"a kg2 list with no bulk prefix should merge on its own, but split into {reps}"


def test_bulk_prefix_ids_stay_out_but_the_rest_of_a_long_list_merges(tmp_path):
    """TP53's shape, end to end: a long kg2 list whose good ids merge while its one bulk prefix does not -- and
    every bulk id is still kept, as its own singleton. DOID/UMLS/NCIT stand in for the good ids; MONDO is avoided
    because its one-id guardrail would split them for the wrong reason."""
    good = ["DOID:1", "UMLS:C1", "NCIT:C1", "MESH:D1", "OMIM:MTHU1", "CHV:1"]
    bulk = [f"REACT:R-HSA-{i}" for i in range(40)]
    config = _aggregator_only_clique_config(tmp_path, [*good, *bulk])  # head = good[0]
    m = resolve_entities(config, biolink=None)
    assert all(m.get(c) is not None for c in good + bulk), "an id was dropped"
    assert len({m[c] for c in good}) == 1, "the list's good ids should merge despite its length"
    for r in bulk:
        assert m[r] == r, f"bulk-prefix id {r} merged; it should only be weakly linked"


def test_every_source_id_survives_as_singleton_with_provenance(tmp_path):
    """A bare equiv-list-only id that Babel doesn't know and that never merges is
    kept as its own singleton node, carrying the referencing source's provenance and
    a category inherited from its referencing node (never silently dropped).

    To get an id that genuinely never merges, WEIRD:1 sits in a bulk prefix (more WEIRD ids than robokop's
    max_ids_per_prefix), so it gets only a weak link. (A two-id list used to serve, back when aggregator lists
    were a flat sub-tau weight; a small list now merges on its own.)"""
    import jsonlines

    padding = [f"WEIRD:{i}" for i in range(2, 13)]  # makes WEIRD a bulk prefix in this list
    node = {
        "id": "MONDO:9",
        "categories": ["biolink:Disease"],
        # WEIRD:1 (bare, unrecognized) is one of 12 WEIRD ids -> bulk, so it won't merge
        "equivalent_ids": ["MONDO:9", "WEIRD:1", *padding],
        "provided_by": ["infores:robokop-kg"],
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"robokop": _write_source(tmp_path, "robokop", [node])},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
    )
    resolve_entities(config, biolink=None)
    by_id = {n["id"]: n for n in jsonlines.open(config.integrated_nodes_path)}
    assert "WEIRD:1" in by_id, "bare equiv-list id was dropped instead of kept as a singleton"
    assert by_id["WEIRD:1"]["provided_by"] == ["infores:robokop-kg"]
    # typed via inheritance from its single-family referencing node (Disease), even
    # though Babel doesn't know it and it never merged
    assert by_id["WEIRD:1"]["categories"] == ["biolink:Disease"]


def test_name_sim_links_ids_by_their_own_babel_labels(tmp_path):
    """Name-similarity operates PER-ID on each id's own Babel label: two ids whose individual labels normalize the
    same are linked, with no equivalence edge and no harmonized primary name involved."""
    babel = _write_source(
        tmp_path,
        "babel",
        [
            _babel_node("DOID:14330", "biolink:Disease", "Parkinson disease"),
            _babel_node("MESH:D010300", "biolink:Disease", "PARKINSON DISEASE"),
        ],
    )
    m = resolve_entities(_config(tmp_path, {"babel": babel}), biolink=None)
    assert m.get("DOID:14330") == m.get("MESH:D010300"), "per-id Babel-label name-sim did not link them"


def test_display_name_and_taxon_from_babel(tmp_path):
    """The merged node's display name is the REPRESENTATIVE id's own Babel label (source name demotes to a
    synonym), and Babel's taxon is retained even though the source node carried none."""
    source_node = {
        "id": "NCBIGene:1636",
        "categories": ["biolink:Gene"],
        "equivalent_ids": ["NCBIGene:1636"],
        "provided_by": ["infores:ncbi-gene"],
        "name": "ACE",  # no taxon on the source node
    }
    harmonized = {
        "ncbigene": _write_source(tmp_path, "ncbigene", [source_node]),
        "babel": _write_source(
            tmp_path,
            "babel",
            [_babel_node("NCBIGene:1636", "biolink:Gene", "Angiotensin I converting enzyme", "NCBITaxon:9606")],
        ),
    }
    config = _config(tmp_path, harmonized)
    resolve_entities(config, biolink=None)
    n = {node["id"]: node for node in jsonlines.open(config.integrated_nodes_path)}["NCBIGene:1636"]
    assert n["name"] == "Angiotensin I converting enzyme"  # rep's Babel label preferred
    assert "ACE" in n.get("synonyms", [])  # source name demoted to synonym
    assert n.get("taxon") == "NCBITaxon:9606"  # Babel taxon retained (source had none)


def test_harmonized_member_babel_label_retained_as_synonym(tmp_path):
    """A HARMONIZED (non-bare) member's own Babel label is retained as a synonym on the merged node -- not just its
    source name (never throw away a label)."""
    ncbigene = [
        {"id": "NCBIGene:1636", "categories": ["biolink:Gene"], "equivalent_ids": ["NCBIGene:1636"],
         "provided_by": ["infores:ncbi-gene"], "name": "ACE"},
        {"id": "HGNC:2707", "categories": ["biolink:Gene"], "equivalent_ids": ["HGNC:2707"],
         "provided_by": ["infores:ncbi-gene"], "name": "ACE_symbol"},
    ]
    babel = [
        _babel_node("NCBIGene:1636", "biolink:Gene", "angiotensin converting enzyme"),
        _babel_node("HGNC:2707", "biolink:Gene", "ACE gene"),
    ]
    harmonized = {
        "ncbigene": _write_source(tmp_path, "ncbigene", ncbigene),
        "babel": _write_source(tmp_path, "babel", babel, _same_as_star(["NCBIGene:1636", "HGNC:2707"])),
    }
    config = _config(tmp_path, harmonized)
    m = resolve_entities(config, biolink=None)
    assert m["NCBIGene:1636"] == m["HGNC:2707"]  # merged via the Babel clique
    merged = {node["id"]: node for node in jsonlines.open(config.integrated_nodes_path)}[m["HGNC:2707"]]
    syns = merged.get("synonyms", [])
    # HGNC:2707 is the rep (name = its Babel label "ACE gene"); NCBIGene:1636's Babel label
    # is a harmonized member's label and must survive as a synonym.
    assert "angiotensin converting enzyme" in syns
    assert "ACE" in syns and "ACE_symbol" in syns  # source names retained too
