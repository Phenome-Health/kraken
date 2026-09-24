"""Rejoining the clusters a Babel gene/protein clique was split across."""

from kraken.entity_resolution.families import ALL_FAMILIES
from kraken.entity_resolution.gene_protein_cohesion import rejoin_split_gene_protein_cliques
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo

GENE_PROTEIN = frozenset({"gene_protein"})
CHEMICAL = frozenset({"chemical"})
HUMAN = "NCBITaxon:9606"

# ACE's shape: one Babel clique, the gene ids in one cluster and the protein ids in another.
CLIQUE = ["HGNC:2707", "NCBIGene:1636", "PR:P12821", "UniProtKB:P12821"]
SPLIT = {"HGNC:2707": 0, "NCBIGene:1636": 0, "PR:P12821": 1, "UniProtKB:P12821": 1}


def _cliques(tmp_path, *cliques):
    path = tmp_path / "cliques.tsv"
    path.write_text("".join(f"{m}\t{c[0]}\t{len(c)}\n" for c in cliques for m in c))
    return path


def _rejoin(tmp_path, clusters, info, cliques=(CLIQUE,)):
    def info_of(curie):
        return info.get(curie, NodeInfo(curie=curie, branches=GENE_PROTEIN, taxon=HUMAN))

    merged = rejoin_split_gene_protein_cliques(clusters, _cliques(tmp_path, *cliques), info_of, GuardrailConfig())
    return clusters, merged


def test_a_gene_and_its_protein_torn_apart_are_put_back(tmp_path):
    clusters, merged = _rejoin(tmp_path, dict(SPLIT), {})
    assert len(set(clusters.values())) == 1
    assert [sorted(cluster) for cluster in merged[0]] == [
        ["HGNC:2707", "NCBIGene:1636"],
        ["PR:P12821", "UniProtKB:P12821"],
    ]


def test_a_clique_already_in_one_cluster_is_left_alone(tmp_path):
    clusters, merged = _rejoin(tmp_path, dict.fromkeys(CLIQUE, 0), {})
    assert merged == []
    assert set(clusters.values()) == {0}


def test_two_species_in_one_clique_stay_apart(tmp_path):
    """The usual refusal: Babel cliques that mix orthologs -- three quarters of the split gene/protein cliques
    in the 2.3.0 build."""
    info = {"PR:P12821": NodeInfo(curie="PR:P12821", branches=GENE_PROTEIN, taxon="NCBITaxon:10090")}
    clusters, merged = _rejoin(tmp_path, dict(SPLIT), info)
    assert merged == []
    assert clusters == SPLIT


def test_a_cluster_holding_anything_but_a_gene_or_protein_is_left_alone(tmp_path):
    """The repair is this family's alone: elsewhere a torn Babel clique is usually label propagation separating
    a class from its members."""
    clusters = {**SPLIT, "CHEBI:15377": 1}
    info = {"CHEBI:15377": NodeInfo(curie="CHEBI:15377", branches=CHEMICAL)}
    clusters, merged = _rejoin(tmp_path, clusters, info)
    assert merged == []
    assert clusters["HGNC:2707"] != clusters["PR:P12821"]


def test_untyped_ids_of_the_clusters_come_along(tmp_path):
    """An untyped id (a LOINC part naming the gene, in ACE's case) is a guardrail wildcard, not a blocker."""
    clusters = {**SPLIT, "LOINC:LP288052-6": 0}
    info = {"LOINC:LP288052-6": NodeInfo(curie="LOINC:LP288052-6", branches=ALL_FAMILIES)}
    clusters, merged = _rejoin(tmp_path, clusters, info)
    assert len(set(clusters.values())) == 1
    assert "LOINC:LP288052-6" in {curie for cluster in merged[0] for curie in cluster}


def test_ids_no_pair_ever_reached_are_not_needed(tmp_path):
    """A clique member with no evidence at all is in no cluster; the rest of its clique still rejoins."""
    clusters, merged = _rejoin(tmp_path, dict(SPLIT), {}, cliques=([*CLIQUE, "UMLS:C1413931"],))
    assert len(set(clusters.values())) == 1
    assert "UMLS:C1413931" not in clusters
