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

from kraken.entity_resolution import build as build_mod
from kraken.entity_resolution.build import resolve_entities
from kraken.entity_resolution.eval.scorer import load_gold, score
from kraken.entity_resolution.sri_nodenorm import NormInfo

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


def _write_source(tmp_path: Path, source: str, records: list[dict]) -> tuple[Path, Path]:
    d = tmp_path / "harmonized" / source
    d.mkdir(parents=True, exist_ok=True)
    nodes_path = d / "nodes.jsonl"
    with jsonlines.open(nodes_path, "w") as w:
        w.write_all(records)
    return nodes_path, d / "edges.jsonl"


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
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
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
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )


class _CliqueNN:
    """Node-normalizer stub that returns the Parkinson clique (the equivalence
    signal now comes from the normalizer's cliques, not aggregator equiv lists)."""

    def __init__(self, *a, **k):
        pass

    def resolve(self, curies, **k):
        return {
            c: NormInfo(label=None, categories=("biolink:Disease",), taxa=(), canonical=PARKINSON[0])
            for c in curies
            if c in PARKINSON
        }

    def iter_cliques(self):
        yield (PARKINSON[0], list(PARKINSON))

    def iter_labels(self):
        return iter(())

    def get(self, curie):
        if curie in PARKINSON:
            return NormInfo(label=None, categories=("biolink:Disease",), canonical=PARKINSON[0])
        return None

    def close(self):
        pass


def test_nn_clique_merges_disease(tmp_path, monkeypatch):
    """A cross-ontology disease clique merges into one node when the NORMALIZER
    supplies the clique (the Parkinson case) -- this is the new equivalence path."""
    monkeypatch.setattr(build_mod, "NodeNormClient", _CliqueNN)
    config = _aggregator_only_clique_config(tmp_path, PARKINSON)
    m = resolve_entities(config, biolink=None)
    mapped = [c for c in PARKINSON if m.get(c) == PARKINSON[0]]
    assert mapped == PARKINSON, (
        f"clique did not fully merge: {len(mapped)}/{len(PARKINSON)} mapped to rep; "
        f"missing/split={[c for c in PARKINSON if m.get(c) != PARKINSON[0]]}"
    )


def test_small_aggregator_equiv_list_merges_on_its_own(tmp_path):
    """A small aggregator list merges without any normalizer clique.

    Aggregator lists are size-aware: small ones are reliable (~95%+ agreement with NN through 20
    ids), so they carry a merge-strength weight. The kg2 Parkinson clique is 10 ids. (Aggregator
    lists used to be a flat sub-tau weight, which left every such list inert -- this test asserted
    the opposite until that was changed.)"""
    config = _aggregator_only_clique_config(tmp_path, PARKINSON)
    m = resolve_entities(config, biolink=None)  # autouse offline NN stub -> no cliques
    reps = {m.get(c) for c in PARKINSON}
    assert None not in reps, "an id was dropped"
    assert len(reps) == 1, f"a 10-id kg2 list should merge on its own, but split into {reps}"


@pytest.mark.parametrize(
    "size, merges",
    [
        (5, True),
        (24, True),  # kg2's max_merge_list_size: still merges
        (25, False),  # one over: corroborates only, below tau
        (60, False),  # max_corroborate_list_size: still corroborates
        (61, False),  # one over: carries no evidence at all
    ],
)
def test_aggregator_equiv_list_merges_on_its_own_only_when_small(tmp_path, size, merges):
    """kg2 merges a list on its own up to 24 ids and not beyond -- and a list that doesn't merge
    still keeps every id, each as its own singleton (nothing a source provided is dropped).

    DOID ids are used deliberately: MONDO is a one-id-per-cluster prefix, so a list of MONDO ids
    would fail to merge because of that guardrail and pass this test for the wrong reason. Here
    size is the only thing that varies."""
    members = [f"DOID:{i}" for i in range(1, size + 1)]
    config = _aggregator_only_clique_config(tmp_path, members)
    m = resolve_entities(config, biolink=None)
    assert all(m.get(c) is not None for c in members), "an id was dropped"
    reps = {m.get(c) for c in members}
    if merges:
        assert len(reps) == 1, f"a {size}-id kg2 list should merge on its own, but split into {len(reps)}"
    else:
        for c in members:
            assert m.get(c) == c, f"a {size}-id kg2 list should NOT merge on its own, but {c} did"


def test_every_source_id_survives_as_singleton_with_provenance(tmp_path):
    """A bare equiv-list-only id that NN doesn't recognize and that never merges is
    kept as its own singleton node, carrying the referencing source's provenance and
    a category inherited from its referencing node (never silently dropped).

    To get an id that genuinely never merges, the robokop list is made one id larger than
    robokop's max_merge_list_size (30), so it only corroborates. (A two-id list used to serve,
    back when aggregator lists were a flat sub-tau weight; a small list now merges on its own.)"""
    import jsonlines

    padding = [f"WEIRDPAD:{i}" for i in range(29)]  # non-enforced prefix, so size alone decides
    node = {
        "id": "MONDO:9",
        "categories": ["biolink:Disease"],
        # 31 ids: past robokop's merge threshold, so WEIRD:1 (bare, unrecognized) won't merge
        "equivalent_ids": ["MONDO:9", "WEIRD:1", *padding],
        "provided_by": ["infores:robokop-kg"],
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"robokop": _write_source(tmp_path, "robokop", [node])},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    resolve_entities(config, biolink=None)
    by_id = {n["id"]: n for n in jsonlines.open(config.integrated_nodes_path)}
    assert "WEIRD:1" in by_id, "bare equiv-list id was dropped instead of kept as a singleton"
    assert by_id["WEIRD:1"]["provided_by"] == ["infores:robokop-kg"]
    # typed via inheritance from its single-family referencing node (Disease), even
    # though NN doesn't recognize it and it never merged
    assert by_id["WEIRD:1"]["categories"] == ["biolink:Disease"]


class _LabelNN:
    """Stub giving two bare ids INDIVIDUAL NN labels that normalize identically."""

    _labels = {"DOID:14330": "Parkinson disease", "MESH:D010300": "PARKINSON DISEASE"}

    def __init__(self, *a, **k):
        pass

    def resolve(self, curies, **k):
        return {
            c: NormInfo(label=self._labels[c], categories=("biolink:Disease",)) for c in curies if c in self._labels
        }

    def iter_cliques(self):
        return iter(())  # no equivalence clique -- ONLY name-sim can link them

    def iter_labels(self):
        return iter(self._labels.items())

    def get(self, curie):
        if curie not in self._labels:
            return None
        return NormInfo(label=self._labels[curie], categories=("biolink:Disease",))

    def close(self):
        pass


def test_name_sim_links_ids_by_individual_nn_label(tmp_path, monkeypatch):
    """Name-similarity operates PER-ID on each id's own NN label: two bare ids whose
    individual NN labels normalize the same are linked, with no equivalence edge and
    no harmonized primary name involved."""
    monkeypatch.setattr(build_mod, "NodeNormClient", _LabelNN)
    node = {
        "id": "DOID:14330",
        "categories": ["biolink:Disease"],
        "equivalent_ids": ["DOID:14330", "MESH:D010300"],
        "provided_by": ["infores:robokop-kg"],
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"robokop": _write_source(tmp_path, "robokop", [node])},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    m = resolve_entities(config, biolink=None)
    assert m.get("DOID:14330") == m.get("MESH:D010300"), "per-id NN-label name-sim did not link them"


class _RepLabelNN:
    _info = NormInfo(
        label="Angiotensin I converting enzyme", categories=("biolink:Gene",), taxa=("NCBITaxon:9606",)
    )

    def __init__(self, *a, **k):
        pass

    def resolve(self, curies, **k):
        return {c: self._info for c in curies if c == "NCBIGene:1636"}

    def iter_cliques(self):
        return iter(())

    def iter_labels(self):
        return iter([("NCBIGene:1636", self._info.label)])

    def get(self, curie):
        return self._info if curie == "NCBIGene:1636" else None

    def close(self):
        pass


def test_display_name_and_taxon_from_normalizer(tmp_path, monkeypatch):
    """The merged node's display name is the REPRESENTATIVE id's own NN label (source
    name demotes to a synonym), and the NN taxon is retained even though the source
    node carried none."""
    monkeypatch.setattr(build_mod, "NodeNormClient", _RepLabelNN)
    node = {
        "id": "NCBIGene:1636",
        "categories": ["biolink:Gene"],
        "equivalent_ids": ["NCBIGene:1636"],
        "provided_by": ["infores:ncbi-gene"],
        "name": "ACE",  # no taxon on the source node
    }
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"ncbigene": _write_source(tmp_path, "ncbigene", [node])},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    resolve_entities(config, biolink=None)
    n = {node["id"]: node for node in jsonlines.open(config.integrated_nodes_path)}["NCBIGene:1636"]
    assert n["name"] == "Angiotensin I converting enzyme"  # rep's NN label preferred
    assert "ACE" in n.get("synonyms", [])  # source name demoted to synonym
    assert n.get("taxon") == "NCBITaxon:9606"  # NN taxon retained (source had none)


class _TwoLabelNN:
    """NN gives each of two harmonized members its own label (differing from source
    name) and groups them in one clique."""

    _labels = {"NCBIGene:1636": "angiotensin converting enzyme", "HGNC:2707": "ACE gene"}

    def __init__(self, *a, **k):
        pass

    def resolve(self, curies, **k):
        return {c: NormInfo(label=self._labels[c], categories=("biolink:Gene",)) for c in curies if c in self._labels}

    def iter_cliques(self):
        yield ("NCBIGene:1636", ["NCBIGene:1636", "HGNC:2707"])

    def iter_labels(self):
        return iter(self._labels.items())

    def get(self, curie):
        return NormInfo(label=self._labels[curie], categories=("biolink:Gene",)) if curie in self._labels else None

    def close(self):
        pass


def test_harmonized_member_nn_label_retained_as_synonym(tmp_path, monkeypatch):
    """A HARMONIZED (non-bare) member's own NN label is retained as a synonym on the
    merged node -- not just its source name (never throw away an NN label)."""
    monkeypatch.setattr(build_mod, "NodeNormClient", _TwoLabelNN)
    nodes = [
        {"id": "NCBIGene:1636", "categories": ["biolink:Gene"], "equivalent_ids": ["NCBIGene:1636"],
         "provided_by": ["infores:ncbi-gene"], "name": "ACE"},
        {"id": "HGNC:2707", "categories": ["biolink:Gene"], "equivalent_ids": ["HGNC:2707"],
         "provided_by": ["infores:ncbi-gene"], "name": "ACE_symbol"},
    ]
    integrated = tmp_path / "integrated"
    config = SimpleNamespace(
        all_harmonized_paths_resolved={"ncbigene": _write_source(tmp_path, "ncbigene", nodes)},
        integrated_dir=integrated,
        integrated_nodes_path=integrated / "nodes.jsonl",
        er_nodenorm_cache_path=tmp_path / "nodenorm.sqlite",
    )
    m = resolve_entities(config, biolink=None)
    assert m["NCBIGene:1636"] == m["HGNC:2707"]  # merged via NN clique
    merged = {node["id"]: node for node in jsonlines.open(config.integrated_nodes_path)}[m["HGNC:2707"]]
    syns = merged.get("synonyms", [])
    # HGNC:2707 is the rep (name = its NN label "ACE gene"); NCBIGene:1636's NN label
    # is a harmonized member's label and must survive as a synonym.
    assert "angiotensin converting enzyme" in syns
    assert "ACE" in syns and "ACE_symbol" in syns  # source names retained too
