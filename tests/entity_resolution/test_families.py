"""Tests for the branch-family map and one-branch logic."""

from pathlib import Path

from kraken.entity_resolution.families import ALL_FAMILIES, BranchFamilies

REPO_ROOT = Path(__file__).resolve().parents[2]
FAMILIES_YAML = REPO_ROOT / "config" / "entity_resolution" / "branch_families.yaml"


def load() -> BranchFamilies:
    return BranchFamilies.load(FAMILIES_YAML)


def test_gene_and_protein_same_family():
    bf = load()
    gene = bf.families_for_category("biolink:Gene")
    protein = bf.families_for_category("biolink:Protein")
    assert gene == protein
    assert "gene_protein" in gene


def test_disease_separate_from_gene_protein():
    bf = load()
    disease = bf.families_for_category("biolink:Disease")
    gene = bf.families_for_category("biolink:Gene")
    assert disease.isdisjoint(gene)


def test_wildcard_categories_are_all():
    bf = load()
    assert bf.families_for_category("biolink:NamedThing") is ALL_FAMILIES
    assert bf.branches(["biolink:NamedThing"]) is ALL_FAMILIES
    assert bf.branches([]) is ALL_FAMILIES
    assert bf.branches(None) is ALL_FAMILIES


def test_clinical_measurement_bridges_but_disease_and_infocontent_stay_apart():
    # The Pompe newborn-screen nodes: LOINC obs = CDE = Procedure/panel must merge,
    # ClinicalFinding bridges to Disease, but generic InformationContentEntity
    # (CommonDataElement) must NOT be mergeable with Disease itself.
    bf = load()
    assert bf.is_valid_cluster([["biolink:CommonDataElement"], ["biolink:ClinicalMeasurement"]])
    assert bf.is_valid_cluster([["biolink:Procedure"], ["biolink:ClinicalMeasurement"]])
    assert bf.is_valid_cluster([["biolink:ClinicalFinding"], ["biolink:Disease"]])  # bridge
    assert not bf.is_valid_cluster([["biolink:Disease"], ["biolink:CommonDataElement"]])  # hard wall
    # a bridge node cannot smuggle the two extremes into one cluster
    assert not bf.is_valid_cluster([["biolink:Disease"], ["biolink:ClinicalFinding"], ["biolink:CommonDataElement"]])


def test_ace_rtd_cluster_is_invalid():
    # The issue #7 conflation: a gene/protein and a disease in one cluster.
    bf = load()
    node_cats = [
        ["biolink:Gene"],  # NCBIGene:1636
        ["biolink:Protein"],  # UniProtKB:P12821
        ["biolink:Disease"],  # MONDO:0017609
    ]
    assert not bf.is_valid_cluster(node_cats)


def test_gene_protein_cluster_is_valid():
    # Gene + Protein together must be allowed (we want gene/protein conflation).
    bf = load()
    assert bf.is_valid_cluster([["biolink:Gene"], ["biolink:Protein"]])


def test_wildcard_node_never_causes_violation():
    bf = load()
    # A disease cluster plus an untyped node stays valid.
    assert bf.is_valid_cluster([["biolink:Disease"], ["biolink:NamedThing"], []])


def test_branches_permissive_within_node_strict_across():
    bf = load()
    # A single node typed both Gene and Disease is permissive within itself...
    within = bf.branches(["biolink:Gene", "biolink:Disease"])
    assert "gene_protein" in within and "disease_pheno" in within
    # ...but across two such single-typed nodes the intersection is empty.
    assert not bf.is_valid_cluster([["biolink:Gene"], ["biolink:Disease"]])
    # The dual-typed node alone is compatible with a pure gene node.
    assert bf.is_valid_cluster([["biolink:Gene", "biolink:Disease"], ["biolink:Gene"]])


def test_cluster_branches_intersection():
    bf = load()
    b_gene = bf.branches(["biolink:Gene"])
    b_disease = bf.branches(["biolink:Disease"])
    assert BranchFamilies.cluster_branches([b_gene, b_disease]) is None
    assert BranchFamilies.cluster_branches([b_gene, b_gene]) == b_gene
    # all-wildcard cluster -> ALL
    assert BranchFamilies.cluster_branches([ALL_FAMILIES, ALL_FAMILIES]) is ALL_FAMILIES
