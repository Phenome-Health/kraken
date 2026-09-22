"""The ER debug database a build writes, and printing a node from it."""

import io

import pytest

from kraken.entity_resolution.build import resolve_entities
from kraken.entity_resolution.debug_db import DebugDb, debug_db_path
from kraken.entity_resolution.inspect_node import render

pytest.importorskip("igraph")


@pytest.fixture
def glycolipid_build(tmp_path):
    from tests.entity_resolution.test_build import _glycolipid_config

    config = _glycolipid_config(tmp_path, ["Glycolipids", "Glycolipids"])
    resolve_entities(config, biolink=None)
    return DebugDb(debug_db_path(config.integrated_dir, None))


def test_every_id_gets_a_row_with_its_node_and_names(glycolipid_build):
    db = glycolipid_build
    node = db.node_of("HMDB:HMDB0302365")
    assert node == db.node_of("CHEBI:33563")
    members = {m.id: m for m in db.members(node)}
    assert {"HMDB:HMDB0302365", "CHEBI:33563", "MESH:D006017", "UMLS:C0017950"} <= set(members)
    hmdb = members["HMDB:HMDB0302365"]
    assert hmdb.babel_name == "Glycolipids"
    assert hmdb.babel_hub == "CHEBI:64583", "its Babel clique is sphingomyelin's, even though it left it"
    assert db.node_of("NOT:THERE") is None


def test_the_evidence_that_moved_it_is_kept(glycolipid_build):
    db = glycolipid_build
    evidence = {(e.a, e.b): e for e in db.evidence_touching(["HMDB:HMDB0302365"])}
    name_match = evidence[("HMDB:HMDB0302365", "MESH:D006017")]
    assert name_match.above_tau and "name_sim" in name_match.kinds
    kg2 = evidence[("CHEBI:33563", "HMDB:HMDB0302365")]
    assert "equiv:kg2" in kg2.kinds, "the restored aggregator claim"


def test_printing_a_node_shows_every_member_and_the_split_clique(glycolipid_build):
    out = io.StringIO()
    render(glycolipid_build, "MESH:D006017", out)
    text = out.getvalue()
    for curie in ("HMDB:HMDB0302365", "CHEBI:33563", "MESH:D006017", "UMLS:C0017950"):
        assert curie in text
    assert "CHEBI:64583" in text and "elsewhere" in text, "the sphingomyelin clique it came from is shown"
    assert "Glycolipids" in text
