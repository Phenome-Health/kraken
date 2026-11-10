"""Tests for canonical-node materialization (order-independent)."""

import random
from pathlib import Path

import pytest

from kraken.entity_resolution.families import BranchFamilies
from kraken.entity_resolution.materialize import PrefixRanking, materialize_cluster

REPO_ROOT = Path(__file__).resolve().parents[2]
CFG = REPO_ROOT / "config" / "entity_resolution"


def _ranking() -> PrefixRanking:
    return PrefixRanking.load(CFG / "prefix_ranking.yaml")


def _families() -> BranchFamilies:
    return BranchFamilies.load(CFG / "branch_families.yaml")


def _gene_members() -> list[dict]:
    return [
        {
            "id": "NCBIGene:1636",
            "categories": ["biolink:Gene"],
            "name": "ACE (from ncbigene)",
            "provided_by": ["ncbigene"],
            "equivalent_ids": ["NCBIGene:1636", "UniProtKB:P12821"],
            "synonyms": ["CD143"],
            "taxon": "NCBITaxon:9606",
        },
        {
            "id": "HGNC:2707",
            "categories": ["biolink:Gene"],
            "name": "ACE",
            "description": "Angiotensin converting enzyme.",
            "provided_by": ["hgnc"],
            "equivalent_ids": ["HGNC:2707"],
            "taxon": "NCBITaxon:9606",
        },
        {
            "id": "UniProtKB:P12821",
            "categories": ["biolink:Protein"],
            "name": "ACE_HUMAN",
            "provided_by": ["uniprot"],
            "equivalent_ids": ["UniProtKB:P12821"],
        },
    ]


def test_representative_is_top_ranked_prefix():
    canonical = materialize_cluster(_gene_members(), _ranking(), _families())
    # gene_protein ranking puts HGNC first
    assert canonical["id"] == "HGNC:2707"
    assert canonical["name"] == "ACE"
    assert canonical["description"] == "Angiotensin converting enzyme."


def test_union_fields():
    canonical = materialize_cluster(_gene_members(), _ranking(), _families())
    assert set(canonical["equivalent_ids"]) >= {"NCBIGene:1636", "HGNC:2707", "UniProtKB:P12821"}
    assert set(canonical["provided_by"]) == {"ncbigene", "hgnc", "uniprot"}
    # other members' names become synonyms
    assert "ACE_HUMAN" in canonical["synonyms"]
    assert "ACE (from ncbigene)" in canonical["synonyms"]
    assert "CD143" in canonical["synonyms"]
    assert canonical["taxon"] == "NCBITaxon:9606"


def test_order_independent():
    members = _gene_members()
    base = materialize_cluster(members, _ranking(), _families())
    for _ in range(20):
        shuffled = members[:]
        random.shuffle(shuffled)
        assert materialize_cluster(shuffled, _ranking(), _families()) == base


def test_description_from_lower_ranked_when_top_lacks_it():
    members = [
        {
            "id": "MONDO:1",
            "categories": ["biolink:Disease"],
            "name": "disease one",
            "provided_by": ["mondo"],
            "equivalent_ids": ["MONDO:1"],
        },  # no description
        {
            "id": "OMIM:2",
            "categories": ["biolink:Disease"],
            "name": "Disease One",
            "description": "A described disease.",
            "provided_by": ["omim"],
            "equivalent_ids": ["OMIM:2"],
        },
    ]
    canonical = materialize_cluster(members, _ranking(), _families())
    assert canonical["id"] == "MONDO:1"  # MONDO ranks above OMIM
    assert canonical["description"] == "A described disease."  # from OMIM


def _member(curie: str, category: str) -> dict:
    return {"id": curie, "categories": [category], "provided_by": ["x"], "equivalent_ids": [curie]}


def test_chemical_holding_an_atc_substance_code_is_a_drug():
    """Tezepelumab: no id of it is Drug-typed, so queries filtered to biolink:Drug missed an approved antibody."""
    members = [_member("UNII:RJ1IW3B4QX", "biolink:ChemicalEntity"), _member("ATC:R03DX11", "biolink:NamedThing")]
    assert "biolink:Drug" in materialize_cluster(members, _ranking(), _families())["categories"]


@pytest.mark.parametrize(
    "members",
    [
        # a shorter code is a drug class
        [_member("CHEBI:1", "biolink:ChemicalEntity"), _member("ATC:R03DX", "biolink:NamedThing")],
        # a few 7-character codes name a group: carbohydrates on the ChEBI class "carbohydrate"
        [_member("CHEBI:16646", "biolink:MolecularEntity"), _member("ATC:B05BA03", "biolink:NamedThing")],
        # the same codes sit on organisms: Lactobacillus is not a drug
        [_member("NCBITaxon:1578", "biolink:OrganismTaxon"), _member("ATC:A07FA01", "biolink:NamedThing")],
        # a code only a member's source list names went to another node
        [{**_member("CHEBI:2", "biolink:ChemicalEntity"), "equivalent_ids": ["CHEBI:2", "ATC:A01AB08"]}],
    ],
    ids=["class code", "group code", "organism", "spilled list"],
)
def test_what_does_not_make_a_drug(members):
    assert "biolink:Drug" not in materialize_cluster(members, _ranking(), _families())["categories"]
