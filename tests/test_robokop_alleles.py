"""ROBOKOP lists an allele's dbSNP rsid as an equivalent id of its CAID node, but an rsid names a POSITION and
each CAID one ALLELE at it (rs7944541 is carried by three CAIDs). So the rsid is moved out of equivalent_ids
and emitted as an `allele member_of position` edge, with one node per position."""

import json

from kraken.harmonizers import base
from kraken.harmonizers.robokop import ALLELE_TO_POSITION_PREDICATE, POSITION_CATEGORY, RobokopHarmonizer
from tests.helpers import PassthroughNormalizer


class _StubBiolink:
    version = "4.2.5"

    def filter_to_leaf_categories(self, categories):
        return list(categories)


def _harmonizer(monkeypatch) -> RobokopHarmonizer:
    """A real RobokopHarmonizer (so __init__ sets up all its state), minus the network: biomapper2's
    Normalizer is swapped for a pass-through stub."""
    monkeypatch.setattr(base, "Normalizer", lambda **_kwargs: PassthroughNormalizer())
    return RobokopHarmonizer(_StubBiolink(), "infores:robokop-kg")


def _write(path, rows):
    path.write_text("".join(json.dumps(r) + "\n" for r in rows))


def _read(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def test_rsid_becomes_a_membership_edge_not_an_equivalent_id(tmp_path, monkeypatch):
    alleles = ["CAID:CA675382683", "CAID:CA1961200538", "CAID:CA220112499"]
    nodes_in = [
        {
            "id": a,
            "name": "rs7944541",
            "category": ["biolink:SequenceVariant"],
            "equivalent_identifiers": [a, "DBSNP:rs7944541", f"HGVS:NC_000011.10:g.{i}A>G"],
        }
        for i, a in enumerate(alleles)
    ] + [{"id": "NCBIGene:1", "name": "g", "category": ["biolink:Gene"], "equivalent_identifiers": ["NCBIGene:1"]}]
    _write(tmp_path / "n.jsonl", nodes_in)
    _write(tmp_path / "e.jsonl", [])

    _harmonizer(monkeypatch).harmonize(
        tmp_path / "nodes.jsonl",
        tmp_path / "edges.jsonl",
        nodes_input=tmp_path / "n.jsonl",
        edges_input=tmp_path / "e.jsonl",
    )
    nodes = {n["id"]: n for n in _read(tmp_path / "nodes.jsonl")}
    edges = _read(tmp_path / "edges.jsonl")

    for a in alleles:
        assert "DBSNP:rs7944541" not in nodes[a]["equivalent_ids"]  # the position is no longer "the allele"
        assert any(e.startswith("HGVS:") for e in nodes[a]["equivalent_ids"])  # but the allele's HGVS still is

    # exactly one position node, however many alleles sit at it
    assert [n for n in nodes if n.startswith("DBSNP:")] == ["DBSNP:rs7944541"]
    assert nodes["DBSNP:rs7944541"]["categories"] == [POSITION_CATEGORY]

    membership = [e for e in edges if e["predicate"] == ALLELE_TO_POSITION_PREDICATE]
    assert sorted(e["subject"] for e in membership) == sorted(alleles)
    assert {e["object"] for e in membership} == {"DBSNP:rs7944541"}

    assert nodes["NCBIGene:1"]["equivalent_ids"] == ["NCBIGene:1"]  # non-allele nodes untouched


def test_smiles_attribute_becomes_a_smiles_equivalent_id(monkeypatch):
    node = _harmonizer(monkeypatch)._harmonize_node(
        {
            "id": "CHEBI:367163",
            "name": "Darunavir",
            "category": ["biolink:SmallMolecule"],
            "equivalent_identifiers": ["CHEBI:367163", "PUBCHEM.COMPOUND:213039"],
            "smiles": "CC(C)CN(C[C@@H](O)[C@H](CC1=CC=CC=C1)NC(=O)O[C@H]1CO[C@H]2OCC[C@@H]12)S(=O)(=O)C1=CC=C(N)C=C1",
        }
    )
    assert set(node["equivalent_ids"]) == {
        "CHEBI:367163",
        "PUBCHEM.COMPOUND:213039",
        "SMILES:CC(C)CN(C[C@@H](O)[C@H](CC1=CC=CC=C1)NC(=O)O[C@H]1CO[C@H]2OCC[C@@H]12)S(=O)(=O)C1=CC=C(N)C=C1",
    }
