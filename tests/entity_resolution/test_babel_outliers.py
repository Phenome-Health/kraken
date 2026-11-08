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


def test_an_id_the_aggregators_place_with_babel_is_left_alone(tmp_path):
    """ChEMBL names CHEMBL455602 "CITPRESSINE II" though its structure is citpressine I, and kg2, ROBOKOP and
    Translator all place it where Babel does -- so the wrong name is ChEMBL's, not Babel's placement."""
    from kraken.entity_resolution.babel_outliers import drop_corroborated

    candidates = {"X": ("H1", "H2", 3)}
    cliques = tmp_path / "cliques.tsv"
    cliques.write_text("".join(f"{m}\tH1\t4\n" for m in ("H1", "A", "B", "X")) + "H2\tH2\t2\nP\tH2\t2\n")
    deferred = tmp_path / "deferred.tsv"

    deferred.write_text("A\tX\tbabel_derived\t0.5\tequiv:kg2\n")  # an aggregator agrees with Babel
    assert drop_corroborated(candidates, cliques, deferred) == {}

    deferred.write_text("P\tX\tbabel_derived\t0.5\tequiv:kg2\n")  # ...but here it says the OTHER clique
    assert drop_corroborated(candidates, cliques, deferred) == candidates

    deferred.write_text("")  # nothing either way: the name still decides
    assert drop_corroborated(candidates, cliques, deferred) == candidates


def test_a_primary_source_mapping_beats_a_list_copied_from_an_older_babel(tmp_path):
    """DOID:1827 "generalized epilepsy" sits in idiopathic generalized epilepsy's clique, kg2's and ROBOKOP's lists copy
    that, and DOID's own xref (a kg2 close_match on the original ids) points at generalized epilepsy's clique."""
    from kraken.entity_resolution.babel_outliers import drop_corroborated

    candidates = {"X": ("H1", "H2", 3)}
    cliques = tmp_path / "cliques.tsv"
    cliques.write_text("".join(f"{m}\tH1\t4\n" for m in ("H1", "A", "B", "X")) + "H2\tH2\t2\nP\tH2\t2\n")
    deferred = tmp_path / "deferred.tsv"
    deferred.write_text("A\tX\tbabel_derived\t0.5\tequiv:kg2\n")
    evidence = tmp_path / "evidence.tsv"

    def kept(rows: str) -> bool:
        evidence.write_text(rows)
        return drop_corroborated(candidates, cliques, deferred, evidence, native_list_sources={"refmet"}) == candidates

    doid = "babel_derived|ks:infores:disease-ontology"
    assert kept(f"P\tX\t{doid}\t0.1\tmatch:kg2\n")

    assert not kept("P\tX\tname_similarity\t0.45\tname_sim\n"), "a name match is what made it a candidate"
    assert not kept("P\tX\tbabel_derived\t0.1\tmatch:kg2\n"), "no primary source named"
    assert not kept("P\tX\tbabel_derived|ks:infores:sri-node-normalizer\t0.1\tmatch:kg2\n"), "Babel confirming itself"
    assert not kept(f"P\tX\t{doid}\t0.05\tmatch_with_subclass:kg2\n"), "kg2 also says subclass: a hierarchy"

    # primary sources disagree: no override, whichever file the edge pointing at its own clique is in
    assert not kept(f"P\tX\t{doid}\t0.1\tmatch:kg2\nA\tX\t{doid}\t0.1\tmatch:kg2\n")
    assert not kept(f"P\tX\t{doid}\t0.1\tmatch:kg2\nB\tX\tsrc:refmet\t1.0\tequiv:refmet\n"), "RefMet's own list"
    assert kept(f"P\tX\t{doid}\t0.1\tmatch:kg2\nB\tX\tbabel_derived\t0.5\tequiv:kg2\n"), "an aggregator's is not"
    assert not kept("P\tX\tsrc:refmet\t1.0\tequiv:refmet\n"), "a native list can keep an id, not move it"
    deferred.write_text(f"A\tX\t{doid}\t1.0\tmatch:kg2\n")
    assert not kept(f"P\tX\t{doid}\t0.1\tmatch:kg2\n")

    # a same_as is set aside with the lists, and counts from there, unless kg2 also says subclass
    deferred.write_text(f"A\tX\tbabel_derived\t0.5\tequiv:kg2\nP\tX\t{doid}\t1.0\tmatch:kg2\n")
    assert kept("")
    deferred.write_text(f"A\tX\tbabel_derived\t0.5\tequiv:kg2\nP\tX\t{doid}\t0.5\tmatch_with_subclass:kg2\n")
    assert not kept("")
