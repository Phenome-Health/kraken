"""Tests for hereditary guardrails and split-until-valid repair."""

from pathlib import Path

from kraken.entity_resolution.families import BranchFamilies
from kraken.entity_resolution.guardrails import (
    DEFAULT_ENFORCED_PREFIXES,
    GuardrailConfig,
    NodeInfo,
    cluster_violations,
    enforce_cluster,
    greedy_valid_partition,
    ids_per_cluster_histogram,
    one_id_valid,
    taxon_valid,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
FAMILIES_YAML = REPO_ROOT / "config" / "entity_resolution" / "branch_families.yaml"


def families() -> BranchFamilies:
    return BranchFamilies.load(FAMILIES_YAML)


def _ace_info() -> dict[str, NodeInfo]:
    return {
        "NCBIGene:1636": NodeInfo("NCBIGene:1636", ("biolink:Gene",), taxon="NCBITaxon:9606"),
        "HGNC:2707": NodeInfo("HGNC:2707", ("biolink:Gene",), taxon="NCBITaxon:9606"),
        "UniProtKB:P12821": NodeInfo("UniProtKB:P12821", ("biolink:Protein",), taxon="NCBITaxon:9606"),
        "MONDO:0017609": NodeInfo("MONDO:0017609", ("biolink:Disease",), taxon=None),
        "orphanet:3033": NodeInfo("orphanet:3033", ("biolink:Disease",), taxon=None),
    }


def test_branch_violation_detected():
    bf, info, cfg = families(), _ace_info(), GuardrailConfig()
    members = list(info)
    assert "branch" in cluster_violations(members, info, bf, cfg)


def test_taxon_guardrail():
    info = {
        "A:1": NodeInfo("A:1", ("biolink:Gene",), taxon="NCBITaxon:9606"),
        "A:2": NodeInfo("A:2", ("biolink:Gene",), taxon="NCBITaxon:10090"),
        "A:3": NodeInfo("A:3", ("biolink:Gene",), taxon=None),  # wildcard taxon
    }
    assert not taxon_valid(["A:1", "A:2"], info)
    assert taxon_valid(["A:1", "A:3"], info)  # untaxoned is a wildcard


def test_one_id_guardrail():
    assert not one_id_valid(["RM:1", "RM:2"], frozenset({"RM"}))
    assert one_id_valid(["RM:1", "CHEBI:2"], frozenset({"RM"}))


def test_default_enforced_prefixes():
    # RefMet=RM, LIPID MAPS=LM (verified in harmonized data) plus MONDO.
    assert DEFAULT_ENFORCED_PREFIXES == frozenset({"RM", "LM", "MONDO"})


def test_default_config_enforces_one_refmet_lipidmaps_mondo():
    bf, cfg = families(), GuardrailConfig()  # defaults enforce RM + LM + MONDO
    info = {
        "RM:1": NodeInfo("RM:1", ("biolink:SmallMolecule",)),
        "RM:2": NodeInfo("RM:2", ("biolink:SmallMolecule",)),
        "LM:1": NodeInfo("LM:1", ("biolink:SmallMolecule",)),
        "LM:2": NodeInfo("LM:2", ("biolink:SmallMolecule",)),
        "MONDO:1": NodeInfo("MONDO:1", ("biolink:Disease",)),
        "MONDO:2": NodeInfo("MONDO:2", ("biolink:Disease",)),
    }
    assert "one_id" in cluster_violations(["RM:1", "RM:2"], info, bf, cfg)
    assert "one_id" in cluster_violations(["LM:1", "LM:2"], info, bf, cfg)
    assert "one_id" in cluster_violations(["MONDO:1", "MONDO:2"], info, bf, cfg)
    # a cluster with two MONDO ids gets split until each part has at most one
    parts = enforce_cluster(["MONDO:1", "MONDO:2"], info, bf, cfg)
    assert len(parts) == 2
    for part in parts:
        assert not cluster_violations(part, info, bf, cfg)


def test_enforce_splits_ace_by_branch():
    bf, info, cfg = families(), _ace_info(), GuardrailConfig()
    # no splitter -> greedy fallback
    parts = enforce_cluster(list(info), info, bf, cfg)
    # every resulting part must be valid
    for part in parts:
        assert not cluster_violations(part, info, bf, cfg)
    # gene/protein land together, disease apart
    gene_part = next(p for p in parts if "NCBIGene:1636" in p)
    assert "UniProtKB:P12821" in gene_part
    assert "MONDO:0017609" not in gene_part


def test_greedy_respects_connectivity_for_wildcards():
    bf, cfg = families(), GuardrailConfig()
    info = {
        "MONDO:1": NodeInfo("MONDO:1", ("biolink:Disease",)),
        "NCBIGene:1": NodeInfo("NCBIGene:1", ("biolink:Gene",)),
        "X:1": NodeInfo("X:1", ("biolink:NamedThing",)),  # wildcard, valid anywhere
    }
    # X:1 is strongly connected to the gene node -> should join the gene group.
    adjacency = {"X:1": {"NCBIGene:1": 5.0, "MONDO:1": 0.1}}
    parts = greedy_valid_partition(list(info), info, bf, cfg, adjacency)
    gene_group = next(p for p in parts if "NCBIGene:1" in p)
    assert "X:1" in gene_group


def test_splitter_used_when_it_reduces():
    bf, info, cfg = families(), _ace_info(), GuardrailConfig()
    calls = {"n": 0}

    def splitter(members):
        calls["n"] += 1
        # pretend Leiden separates gene/protein from disease
        genes = [m for m in members if m.split(":")[0] in {"NCBIGene", "HGNC", "UniProtKB"}]
        disease = [m for m in members if m.split(":")[0] in {"MONDO", "orphanet"}]
        return [genes, disease] if genes and disease else [members]

    parts = enforce_cluster(list(info), info, bf, cfg, splitter=splitter)
    assert calls["n"] >= 1
    for part in parts:
        assert not cluster_violations(part, info, bf, cfg)


def test_one_id_repair_cap_leaves_intact():
    cfg = GuardrailConfig(enforced_prefixes=frozenset({"HGNC"}), one_id_repair_cap=3)
    bf = families()
    info = {f"HGNC:{i}": NodeInfo(f"HGNC:{i}", ("biolink:Gene",)) for i in range(6)}
    # 6 HGNC ids > cap 3 -> intact + logged
    parts = enforce_cluster(list(info), info, bf, cfg)
    assert len(parts) == 1


def test_histogram():
    clusters = [["HGNC:1", "HGNC:2", "NCBIGene:1"], ["HGNC:3"]]
    hist = ids_per_cluster_histogram(clusters)
    assert hist["HGNC"] == {2: 1, 1: 1}
    assert hist["NCBIGene"] == {1: 1}
