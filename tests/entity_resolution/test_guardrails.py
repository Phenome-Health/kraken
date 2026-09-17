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


def _complete(members) -> dict[str, dict[str, float]]:
    """Adjacency linking every pair of members at weight 1.0 -- a cluster whose ids all vouch for each other."""
    return {a: {b: 1.0 for b in members if b != a} for a in members}


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
    # every allele linked to the shared rsid, which is what pulled them together
    adjacency = {"DBSNP:rs7944541": {c: 1.0 for c in info if c.startswith("CAID:")}}
    for c in info:
        if c.startswith("CAID:"):
            adjacency[c] = {"DBSNP:rs7944541": 1.0}
    parts = enforce_cluster(list(info), info, cfg, adjacency=adjacency)
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
    parts = enforce_cluster(["MONDO:1", "MONDO:2"], info, cfg, adjacency=_complete(["MONDO:1", "MONDO:2"]))
    assert len(parts) == 2
    for part in parts:
        assert not cluster_violations(part, info, cfg)


def test_enforce_splits_ace_by_branch():
    info, cfg = _ace_info(), GuardrailConfig()
    # the gene/protein ids vouch for each other, the disease ids for each other, and one conflation edge bridges them
    genes = ["NCBIGene:1636", "HGNC:2707", "UniProtKB:P12821"]
    adjacency = {**_complete(genes), **{d: {} for d in ("MONDO:0017609", "orphanet:3033")}}
    adjacency["MONDO:0017609"]["orphanet:3033"] = adjacency["orphanet:3033"]["MONDO:0017609"] = 1.0
    adjacency["MONDO:0017609"]["NCBIGene:1636"] = adjacency["NCBIGene:1636"]["MONDO:0017609"] = 0.3
    parts = enforce_cluster(list(info), info, cfg, adjacency=adjacency)  # no splitter -> greedy fallback
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


def test_split_keeps_the_strongest_link_rather_than_following_id_order():
    """MONDO:1 -0.5- UMLS:C1 -0.5- MONDO:2 -1.0- DOID:7 holds two MONDO ids, so it must split. Placing ids in sorted
    order put MONDO:1 in DOID:7's group before MONDO:2 (DOID:7's only real link) was placed, cutting that link."""
    cfg = GuardrailConfig()
    members = ["MONDO:1", "UMLS:C1", "MONDO:2", "DOID:7"]
    info = {m: _ni(m, ("biolink:Disease",)) for m in members}
    adjacency: dict[str, dict[str, float]] = {m: {} for m in members}
    for a, b, weight in [("MONDO:1", "UMLS:C1", 0.5), ("UMLS:C1", "MONDO:2", 0.5), ("MONDO:2", "DOID:7", 1.0)]:
        adjacency[a][b] = adjacency[b][a] = weight
    parts = greedy_valid_partition(members, info, cfg, adjacency)
    assert ["DOID:7", "MONDO:2"] in parts
    assert all(not cluster_violations(part, info, cfg) for part in parts)


def test_split_never_groups_ids_without_an_edge_between_them():
    """A member whose every edge would break a guardrail stays on its own -- it isn't parked in some group it has no
    evidence for."""
    cfg = GuardrailConfig()
    members = ["MONDO:1", "MONDO:2", "DOID:9"]
    info = {m: _ni(m, ("biolink:Disease",)) for m in members}
    adjacency = {"MONDO:1": {"MONDO:2": 1.0}, "MONDO:2": {"MONDO:1": 1.0}, "DOID:9": {}}
    parts = greedy_valid_partition(members, info, cfg, adjacency)
    assert sorted(parts) == [["DOID:9"], ["MONDO:1"], ["MONDO:2"]]


def test_splitter_used_when_it_reduces():
    info, cfg = _ace_info(), GuardrailConfig()
    calls = {"n": 0}

    def splitter(members):
        calls["n"] += 1
        genes = [m for m in members if m.split(":")[0] in {"NCBIGene", "HGNC", "UniProtKB"}]
        disease = [m for m in members if m.split(":")[0] in {"MONDO", "orphanet"}]
        return [genes, disease] if genes and disease else [members]

    parts = enforce_cluster(list(info), info, cfg, adjacency=_complete(list(info)), splitter=splitter)
    assert calls["n"] >= 1
    for part in parts:
        assert not cluster_violations(part, info, cfg)


def test_large_one_id_violation_is_repaired_and_reported(caplog):
    """k ids of a one-entity-per-id prefix means k merged entities, however large k is.

    This used to be capped: past 3 ids the cluster was left intact, which shipped the worst conflations
    (a 2.1.1 cluster with 247 RefMet ids) unrepaired. Large repairs are now carried out, and tallied for the
    caller to report -- per-cluster warnings drowned the log, since a build does millions of these."""
    from collections import Counter

    from kraken.entity_resolution.guardrails import log_one_id_repairs

    cfg = GuardrailConfig(enforced_prefixes=frozenset({"HGNC"}), one_id_repair_log_threshold=3)
    info = {f"HGNC:{i}": _ni(f"HGNC:{i}", ("biolink:Gene",)) for i in range(6)}
    repairs: Counter = Counter()
    parts = enforce_cluster(list(info), info, cfg, adjacency=_complete(list(info)), repairs=repairs)
    assert len(parts) == 6
    assert all(not cluster_violations(part, info, cfg) for part in parts)
    assert repairs["HGNC"] == 1  # one cluster repaired, however many ids it held

    with caplog.at_level("INFO"):
        log_one_id_repairs(repairs, cfg)
    assert "split 1 clusters holding more than one HGNC id" in caplog.text
    assert "worst: 6" in caplog.text


def test_a_small_one_id_repair_is_counted_but_not_called_out(caplog):
    from collections import Counter

    from kraken.entity_resolution.guardrails import log_one_id_repairs

    cfg = GuardrailConfig(enforced_prefixes=frozenset({"HGNC"}), one_id_repair_log_threshold=3)
    info = {f"HGNC:{i}": _ni(f"HGNC:{i}", ("biolink:Gene",)) for i in range(2)}
    repairs: Counter = Counter()
    parts = enforce_cluster(list(info), info, cfg, adjacency=_complete(list(info)), repairs=repairs)
    assert len(parts) == 2
    assert repairs["HGNC"] == 1
    with caplog.at_level("INFO"):
        log_one_id_repairs(repairs, cfg)
    assert "upstream conflation" not in caplog.text  # only repairs over the threshold are called out


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
