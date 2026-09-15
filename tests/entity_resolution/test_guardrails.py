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
BF = BranchFamilies.load(REPO_ROOT / "config" / "entity_resolution" / "branch_families.yaml")


def _ni(curie: str, categories: tuple[str, ...] = (), taxon: str | None = None) -> NodeInfo:
    """NodeInfo carries the resolved branch-family set; build it from categories."""
    return NodeInfo(curie=curie, branches=BF.branches(categories), taxon=taxon)


def _ace_info() -> dict[str, NodeInfo]:
    return {
        "NCBIGene:1636": _ni("NCBIGene:1636", ("biolink:Gene",), taxon="NCBITaxon:9606"),
        "HGNC:2707": _ni("HGNC:2707", ("biolink:Gene",), taxon="NCBITaxon:9606"),
        "UniProtKB:P12821": _ni("UniProtKB:P12821", ("biolink:Protein",), taxon="NCBITaxon:9606"),
        "MONDO:0017609": _ni("MONDO:0017609", ("biolink:Disease",), taxon=None),
        "orphanet:3033": _ni("orphanet:3033", ("biolink:Disease",), taxon=None),
    }


def test_branch_violation_detected():
    info, cfg = _ace_info(), GuardrailConfig()
    assert "branch" in cluster_violations(list(info), info, cfg)


def test_taxon_guardrail():
    info = {
        "A:1": _ni("A:1", ("biolink:Gene",), taxon="NCBITaxon:9606"),
        "A:2": _ni("A:2", ("biolink:Gene",), taxon="NCBITaxon:10090"),
        "A:3": _ni("A:3", ("biolink:Gene",), taxon=None),  # wildcard taxon
    }
    assert not taxon_valid(["A:1", "A:2"], info)
    assert taxon_valid(["A:1", "A:3"], info)  # untaxoned is a wildcard


def test_one_id_guardrail():
    assert not one_id_valid(["RM:1", "RM:2"], frozenset({"RM"}))
    assert one_id_valid(["RM:1", "CHEBI:2"], frozenset({"RM"}))


def test_default_enforced_prefixes():
    # RefMet=RM, LIPID MAPS=LM (verified in harmonized data), MONDO, ClinGen alleles (CAID), and structures (SMILES).
    assert DEFAULT_ENFORCED_PREFIXES == frozenset({"RM", "LM", "MONDO", "CAID", "SMILES"})


def test_two_alleles_are_never_one_cluster():
    """A CAID is one allele. Two in a cluster means distinct alleles were merged -- e.g. through a shared
    rsid, which names the position they sit at, not either allele."""
    cfg = GuardrailConfig()
    info = {
        c: _ni(c, ("biolink:SequenceVariant",))
        for c in ("CAID:CA675382683", "CAID:CA1961200538", "CAID:CA220112499", "DBSNP:rs7944541")
    }
    parts = enforce_cluster(list(info), info, cfg)
    assert all(sum(m.startswith("CAID:") for m in part) <= 1 for part in parts)


def test_default_config_enforces_one_refmet_lipidmaps_mondo():
    cfg = GuardrailConfig()  # defaults enforce RM + LM + MONDO
    info = {
        "RM:1": _ni("RM:1", ("biolink:SmallMolecule",)),
        "RM:2": _ni("RM:2", ("biolink:SmallMolecule",)),
        "LM:1": _ni("LM:1", ("biolink:SmallMolecule",)),
        "LM:2": _ni("LM:2", ("biolink:SmallMolecule",)),
        "MONDO:1": _ni("MONDO:1", ("biolink:Disease",)),
        "MONDO:2": _ni("MONDO:2", ("biolink:Disease",)),
    }
    assert "one_id" in cluster_violations(["RM:1", "RM:2"], info, cfg)
    assert "one_id" in cluster_violations(["LM:1", "LM:2"], info, cfg)
    assert "one_id" in cluster_violations(["MONDO:1", "MONDO:2"], info, cfg)
    # a cluster with two MONDO ids gets split until each part has at most one
    parts = enforce_cluster(["MONDO:1", "MONDO:2"], info, cfg)
    assert len(parts) == 2
    for part in parts:
        assert not cluster_violations(part, info, cfg)


def test_enforce_splits_ace_by_branch():
    info, cfg = _ace_info(), GuardrailConfig()
    parts = enforce_cluster(list(info), info, cfg)  # no splitter -> greedy fallback
    for part in parts:
        assert not cluster_violations(part, info, cfg)  # every resulting part valid
    # gene/protein land together, disease apart
    gene_part = next(p for p in parts if "NCBIGene:1636" in p)
    assert "UniProtKB:P12821" in gene_part
    assert "MONDO:0017609" not in gene_part


def test_greedy_respects_connectivity_for_wildcards():
    cfg = GuardrailConfig()
    info = {
        "MONDO:1": _ni("MONDO:1", ("biolink:Disease",)),
        "NCBIGene:1": _ni("NCBIGene:1", ("biolink:Gene",)),
        "X:1": _ni("X:1", ("biolink:NamedThing",)),  # wildcard, valid anywhere
    }
    # X:1 is strongly connected to the gene node -> should join the gene group.
    adjacency = {"X:1": {"NCBIGene:1": 5.0, "MONDO:1": 0.1}}
    parts = greedy_valid_partition(list(info), info, cfg, adjacency)
    gene_group = next(p for p in parts if "NCBIGene:1" in p)
    assert "X:1" in gene_group


def test_splitter_used_when_it_reduces():
    info, cfg = _ace_info(), GuardrailConfig()
    calls = {"n": 0}

    def splitter(members):
        calls["n"] += 1
        genes = [m for m in members if m.split(":")[0] in {"NCBIGene", "HGNC", "UniProtKB"}]
        disease = [m for m in members if m.split(":")[0] in {"MONDO", "orphanet"}]
        return [genes, disease] if genes and disease else [members]

    parts = enforce_cluster(list(info), info, cfg, splitter=splitter)
    assert calls["n"] >= 1
    for part in parts:
        assert not cluster_violations(part, info, cfg)


def test_large_one_id_violation_is_repaired_and_logged(caplog):
    """k ids of a one-entity-per-id prefix means k merged entities, however large k is.

    This used to be capped: past 3 ids the cluster was left intact, which shipped the worst conflations
    (a 2.1.1 cluster with 247 RefMet ids) unrepaired. Large repairs are now carried out, and logged."""
    cfg = GuardrailConfig(enforced_prefixes=frozenset({"HGNC"}), one_id_repair_log_threshold=3)
    info = {f"HGNC:{i}": _ni(f"HGNC:{i}", ("biolink:Gene",)) for i in range(6)}
    with caplog.at_level("WARNING"):
        parts = enforce_cluster(list(info), info, cfg)
    assert len(parts) == 6
    assert all(not cluster_violations(part, info, cfg) for part in parts)
    assert "one_id violation with 6 ids" in caplog.text


def test_small_one_id_repair_is_not_logged(caplog):
    cfg = GuardrailConfig(enforced_prefixes=frozenset({"HGNC"}), one_id_repair_log_threshold=3)
    info = {f"HGNC:{i}": _ni(f"HGNC:{i}", ("biolink:Gene",)) for i in range(2)}
    with caplog.at_level("WARNING"):
        parts = enforce_cluster(list(info), info, cfg)
    assert len(parts) == 2
    assert "one_id violation" not in caplog.text


def test_histogram():
    clusters = [["HGNC:1", "HGNC:2", "NCBIGene:1"], ["HGNC:3"]]
    hist = ids_per_cluster_histogram(clusters)
    assert hist["HGNC"] == {2: 1, 1: 1}
    assert hist["NCBIGene"] == {1: 1}


def test_a_compound_and_its_salt_are_never_one_cluster():
    """Canonical SMILES name one structure each; the same lipid from lipidmaps and its sodium salt from translator
    (lumped under the parent's InChIKey) must not share a cluster."""
    assert not one_id_valid(["SMILES:O=C(O)CCCO", "SMILES:O=C([O-])CCCO.[Na+]"], DEFAULT_ENFORCED_PREFIXES)
    assert one_id_valid(["SMILES:O=C(O)CCCO", "LM:FA01050006"], DEFAULT_ENFORCED_PREFIXES)
