"""Finding Babel clique outliers from name matches and clique membership."""

from kraken.entity_resolution.babel_outliers import find_babel_outliers

CHEMICAL = "chemical"


def _outliers(tmp_path, name_pairs, cliques, clique_cap=100):
    pairs_path = tmp_path / "pairs.tsv"
    pairs_path.write_text("".join(f"{a}\t{b}\t{fa}\t{fb}\n" for a, b, fa, fb in name_pairs))
    cliques_path = tmp_path / "cliques.tsv"
    sizes = {hub: len(members) for hub, members in cliques.items()}
    cliques_path.write_text("".join(f"{m}\t{hub}\t{sizes[hub]}\n" for hub, members in cliques.items() for m in members))
    return find_babel_outliers(pairs_path, cliques_path, clique_cap, tmp_path)


# X sits in clique H1, whose other members A and B share a name; X's name matches P and Q of clique H2.
CLIQUES = {"H1": ["H1", "A", "B", "X"], "H2": ["H2", "P", "Q"]}
COHERENT = [("A", "B", CHEMICAL, CHEMICAL)]
MATCHES_ELSEWHERE = [("P", "X", CHEMICAL, CHEMICAL), ("Q", "X", CHEMICAL, CHEMICAL)]


def test_an_id_named_like_another_clique_and_unlike_its_own_is_an_outlier(tmp_path):
    assert _outliers(tmp_path, COHERENT + MATCHES_ELSEWHERE, CLIQUES) == {"X": ("H1", "H2", 2)}


def test_not_an_outlier_if_its_name_matches_its_own_clique(tmp_path):
    pairs = COHERENT + MATCHES_ELSEWHERE + [("A", "X", CHEMICAL, CHEMICAL)]
    assert _outliers(tmp_path, pairs, CLIQUES) == {}


def test_not_an_outlier_if_its_own_clique_agrees_on_no_name(tmp_path):
    """Without a shared name inside the clique there's nothing to be an outlier from: its members may just be
    named differently (a systematic name vs a trivial one)."""
    assert _outliers(tmp_path, MATCHES_ELSEWHERE, CLIQUES) == {}


def test_not_an_outlier_on_a_single_match_elsewhere(tmp_path):
    assert _outliers(tmp_path, COHERENT + MATCHES_ELSEWHERE[:1], CLIQUES) == {}


def test_genes_organisms_and_oversized_cliques_are_not_judged(tmp_path):
    gene = [(a, b, "gene_protein", "gene_protein") for a, b, *_ in COHERENT + MATCHES_ELSEWHERE]
    assert _outliers(tmp_path, gene, CLIQUES) == {}
    # a clique over the cap is a star from its hub: dropping one id could disconnect it
    assert _outliers(tmp_path, COHERENT + MATCHES_ELSEWHERE, CLIQUES, clique_cap=3) == {}
