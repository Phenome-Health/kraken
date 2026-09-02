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

pytest.importorskip("leidenalg")

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
