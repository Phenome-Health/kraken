"""RefMet maps individual lipid species to the KEGG entry for their whole class.

843 glucosylceramide species all carry KEGG.COMPOUND:C01190 ("Glucosylceramide"). Listed as an equivalent id,
that claims each species IS the class, and in 2.1.1 it fused hundreds of species into one node. An xref shared
by more than MAX_ENTRIES_PER_XREF entries is therefore moved out of equivalent_ids into an attribute.
"""

from kraken.harmonizers.refmet import CLASS_XREFS_ATTRIBUTE, MAX_ENTRIES_PER_XREF, RefMetHarmonizer


def _harmonizer() -> RefMetHarmonizer:
    h = object.__new__(RefMetHarmonizer)  # skip __init__ (Biolink toolkit + Normalizer reach the network)
    h.source_infores = "infores:refmet"
    return h


def _node(rm_id, *xrefs):
    return {"id": rm_id, "name": rm_id, "equivalent_ids": [rm_id, *xrefs], "attributes": {"infores:refmet": {}}}


def _nodes(*nodes):
    return {n["id"]: n for n in nodes}


def test_xref_shared_by_many_entries_is_moved_to_an_attribute():
    species = [
        _node(f"RM:{i}", "KEGG.COMPOUND:C01190", f"PUBCHEM.COMPOUND:{i}") for i in range(MAX_ENTRIES_PER_XREF + 1)
    ]
    nodes = _nodes(*species)
    _harmonizer()._strip_class_level_xrefs(nodes)
    for i, node in enumerate(nodes.values()):
        assert node["equivalent_ids"] == [f"RM:{i}", f"PUBCHEM.COMPOUND:{i}"]  # class id gone, the rest intact
        assert node["attributes"]["infores:refmet"][CLASS_XREFS_ATTRIBUTE] == ["KEGG.COMPOUND:C01190"]


def test_xref_shared_by_up_to_the_limit_is_left_alone():
    """Small overlaps (RefMet duplicates, stereo variants) are for entity resolution, not this filter."""
    species = [_node(f"RM:{i}", "CHEBI:1") for i in range(MAX_ENTRIES_PER_XREF)]
    nodes = _nodes(*species)
    _harmonizer()._strip_class_level_xrefs(nodes)
    for node in nodes.values():
        assert "CHEBI:1" in node["equivalent_ids"]
        assert CLASS_XREFS_ATTRIBUTE not in node["attributes"]["infores:refmet"]


def test_one_to_one_xrefs_are_untouched():
    nodes = _nodes(_node("RM:1", "KEGG.COMPOUND:C10565"), _node("RM:2", "KEGG.COMPOUND:C09020"))
    _harmonizer()._strip_class_level_xrefs(nodes)
    assert nodes["RM:1"]["equivalent_ids"] == ["RM:1", "KEGG.COMPOUND:C10565"]
    assert nodes["RM:2"]["equivalent_ids"] == ["RM:2", "KEGG.COMPOUND:C09020"]
