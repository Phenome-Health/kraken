"""Tests for aggregator edge un-canonicalization (shared by match evidence + edge remap)."""

from kraken.entity_resolution.uncanonicalize import original_endpoints
from kraken.utils.constants import EDGE_ATTRIBUTES, EDGE_OBJECT, EDGE_SUBJECT

# KG2's per-source attribute key is its build_config source_id; the value used here is just a realistic
# stand-in, since un-canonicalization now finds kg2pre_ids regardless of which infores key it sits under.
KG2_INFORES = "infores:rtx-kg2"


def test_native_source_uses_own_endpoints():
    edge = {EDGE_SUBJECT: "RM:1", EDGE_OBJECT: "CHEBI:2"}
    assert original_endpoints(edge, "refmet") == [("RM:1", "CHEBI:2")]


def test_kg2_recovers_multiple_original_pairs():
    # kg2 stores canonical endpoints; the originals (possibly several) live in kg2pre_ids
    # as "subject---rel---q---q---q---object---src". One merged edge -> several edges.
    edge = {
        EDGE_SUBJECT: "NCBIGene:7157",
        EDGE_OBJECT: "MONDO:1",
        EDGE_ATTRIBUTES: {
            KG2_INFORES: {
                "kg2pre_ids": [
                    "UniProtKB:P04637---affects---q---q---q---MONDO:1---kg2",
                    "HGNC:11998---affects---q---q---q---DOID:2---kg2",
                ]
            }
        },
    }
    pairs = original_endpoints(edge, "kg2")
    assert pairs == [("UniProtKB:P04637", "MONDO:1"), ("HGNC:11998", "DOID:2")]


def test_robokop_uses_original_subject_object_attribute():
    # robokop/translator/microbiome/multiomics store a single original pair in their
    # per-source attribute dict (under whatever infores key), which we recover.
    edge = {
        EDGE_SUBJECT: "NCBIGene:99",  # Babel-canonical
        EDGE_OBJECT: "MONDO:1",
        EDGE_ATTRIBUTES: {"infores:robokop-kg": {"original_subject": "CAID:CA1", "original_object": "MONDO:1"}},
    }
    assert original_endpoints(edge, "robokop") == [("CAID:CA1", "MONDO:1")]


def test_translator_uses_original_attribute_under_its_own_key():
    edge = {
        EDGE_SUBJECT: "NCBIGene:98558",
        EDGE_OBJECT: "UBERON:1",
        EDGE_ATTRIBUTES: {"translator-kg-open": {"original_subject": "ENSEMBL:X", "original_object": "UBERON:1"}},
    }
    assert original_endpoints(edge, "translator-kg-open") == [("ENSEMBL:X", "UBERON:1")]


def test_canonicalized_aggregator_without_originals_returns_none():
    # no kg2pre_ids and no original_subject/object -> None (caller falls back to stored).
    assert original_endpoints({EDGE_SUBJECT: "A:1", EDGE_OBJECT: "B:2"}, "robokop") is None
    assert original_endpoints({EDGE_SUBJECT: "A:1", EDGE_OBJECT: "B:2"}, "kg2") is None


def test_kg2_falls_back_to_original_attribute_without_pre_ids():
    edge = {
        EDGE_SUBJECT: "NCBIGene:7157",
        EDGE_OBJECT: "MONDO:1",
        EDGE_ATTRIBUTES: {KG2_INFORES: {"original_subject": "UniProtKB:P04637", "original_object": "MONDO:1"}},
    }
    assert original_endpoints(edge, "kg2") == [("UniProtKB:P04637", "MONDO:1")]


def test_original_alias_pairs_maps_each_original_to_its_stored_endpoint():
    """An aggregator edge asserts 'I resolved X to Y' for each endpoint."""
    from kraken.entity_resolution.uncanonicalize import original_alias_pairs

    edge = {
        "subject": "CAID:CA1",
        "object": "NCBIGene:5",
        "attributes": {"infores:robokop-kg": {"original_subject": "HGVS:NC_1:g.1A>G", "original_object": "ENSEMBL:E1"}},
    }
    assert original_alias_pairs(edge, "robokop") == [
        ("HGVS:NC_1:g.1A>G", "CAID:CA1"),
        ("ENSEMBL:E1", "NCBIGene:5"),
    ]


def test_original_alias_pairs_skips_originals_equal_to_the_stored_id():
    """Nothing is asserted when the aggregator did not actually rewrite the endpoint."""
    from kraken.entity_resolution.uncanonicalize import original_alias_pairs

    edge = {
        "subject": "CAID:CA1",
        "object": "NCBIGene:5",
        "attributes": {"infores:robokop-kg": {"original_subject": "CAID:CA1", "original_object": "ENSEMBL:E1"}},
    }
    assert original_alias_pairs(edge, "robokop") == [("ENSEMBL:E1", "NCBIGene:5")]


def test_original_alias_pairs_covers_every_kg2_pre_id_pair():
    """One merged KG2 edge can carry several originals; each contributes both aliases."""
    from kraken.entity_resolution.uncanonicalize import original_alias_pairs

    edge = {
        "subject": "UNII:1",
        "object": "PUBCHEM.COMPOUND:1",
        "attributes": {
            "infores:rtx-kg2": {
                "kg2pre_ids": [
                    "ATC:X---rel---None---None---None---UMLS:Y---src",
                    "CHEBI:Z---rel---None---None---None---UMLS:Y---src",
                ]
            }
        },
    }
    assert original_alias_pairs(edge, "kg2") == [
        ("ATC:X", "UNII:1"),
        ("UMLS:Y", "PUBCHEM.COMPOUND:1"),
        ("CHEBI:Z", "UNII:1"),
        ("UMLS:Y", "PUBCHEM.COMPOUND:1"),
    ]


def test_original_alias_pairs_requires_both_stored_endpoints():
    from kraken.entity_resolution.uncanonicalize import original_alias_pairs

    edge = {"subject": "CAID:CA1", "attributes": {"x": {"original_subject": "HGVS:1", "original_object": "E:1"}}}
    assert original_alias_pairs(edge, "robokop") == []
