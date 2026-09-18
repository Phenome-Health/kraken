"""Metagraph node statistics: taxon counts, named from the graph's own OrganismTaxon nodes."""

import json

from kraken.metagraph import generate_metagraph_streaming


def _write(path, rows):
    path.write_text("".join(json.dumps(row) + "\n" for row in rows))


def test_node_taxa_are_counted_and_named_from_organism_taxon_nodes(tmp_path):
    def node(curie, category, taxon=None, name=None, equivalent_ids=None):
        row = {"id": curie, "categories": [category], "equivalent_ids": equivalent_ids or [curie]}
        if taxon:
            row["taxon"] = taxon
        if name:
            row["name"] = name
        return row

    _write(
        tmp_path / "nodes.jsonl",
        [
            node("HGNC:2707", "biolink:Gene", "NCBITaxon:9606"),
            node("UniProtKB:P12821", "biolink:Protein", "NCBITaxon:9606"),
            node("RGD:2493", "biolink:Gene", "NCBITaxon:10116"),
            node("MONDO:0005148", "biolink:Disease"),
            # the human taxon node's representative isn't the NCBITaxon id -- names still resolve via equivalent_ids
            node(
                "UMLS:C0086418",
                "biolink:OrganismTaxon",
                name="Homo sapiens",
                equivalent_ids=["UMLS:C0086418", "NCBITaxon:9606"],
            ),
        ],
    )
    _write(tmp_path / "edges.jsonl", [])
    stats = generate_metagraph_streaming(tmp_path / "nodes.jsonl", tmp_path / "edges.jsonl", "kraken", "test")
    metagraph = stats.to_dict()

    assert metagraph["node_taxa"] == {
        "NCBITaxon:9606": {"name": "Homo sapiens", "count": 2},
        "NCBITaxon:10116": {"name": None, "count": 1},  # no rat taxon node in this graph
    }
    assert metagraph["summary"]["unique_node_taxa"] == 2
    assert metagraph["summary"]["nodes_without_taxon"] == 2  # the disease and the taxon node itself
