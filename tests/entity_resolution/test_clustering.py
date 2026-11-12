"""Tests for the clustering step: the clusters are a function of the evidence only (see clustering)."""

import json
import random
from collections import defaultdict
from fractions import Fraction
from functools import reduce
from itertools import combinations

from kraken.entity_resolution.clustering import Tie, cluster_component, log_ties
from kraken.entity_resolution.families import BranchFamilies
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo, cluster_violations, group_state, join_states
from kraken.entity_resolution.weights import WEIGHT_UNIT, weight_units

CONFIG = GuardrailConfig()
FAMILIES = BranchFamilies.load()
TAU = weight_units(0.3)
BABEL_WEIGHT = 0.5
HUMAN, MOUSE = "NCBITaxon:9606", "NCBITaxon:10090"
# Any two of these share a family (exposure, chemical, gene_protein), and the three share none.
FAMILY_TRIANGLE = ("biolink:ChemicalExposure", "biolink:GenomicBackgroundExposure", "biolink:NucleicAcidEntity")


def _info(curie, category=None, taxon=None):
    return NodeInfo(curie=curie, branches=FAMILIES.branches((category,)), taxon=taxon)


def _cluster(pairs, cliques=(), taxa=None, categories=None):
    """Cluster what stage 3 would hand over for these pairs and Babel cliques: each clique's own pairs added, and
    the pairs a guardrail forbids pruned."""
    pairs = [*pairs, *((a, b, BABEL_WEIGHT) for clique in cliques for a in clique for b in clique if a < b)]
    members = sorted({m for a, b, _weight in pairs for m in (a, b)})
    info = {m: _info(m, (categories or {}).get(m), (taxa or {}).get(m)) for m in members}
    kept = [(a, b, weight_units(weight)) for a, b, weight in pairs if not cluster_violations([a, b], info, CONFIG)]
    clique_of = {m: number for number, clique in enumerate(cliques) for m in clique}
    return cluster_component(members, kept, info, CONFIG, clique_of, TAU)


def _random_component(rng, label="", first_clique=0):
    """What ``cluster_component`` takes, at random: ids with random types, taxa and enforced prefixes, Babel cliques,
    and pairs on a coarse weight grid, half the time with one id linked equally to several others, so that exact
    ties, chains of equal links and guardrail conflicts are common."""
    members = [f"{rng.choice(['MONDO', 'CHEBI', 'UMLS', 'MESH'])}:{label}{i}" for i in range(rng.randint(2, 12))]
    categories = [None, None, None, "biolink:SmallMolecule", *FAMILY_TRIANGLE]
    info = {m: _info(m, rng.choice(categories), rng.choice([None, None, None, HUMAN, MOUSE])) for m in members}
    shuffled = rng.sample(members, len(members))
    cliques = [shuffled[i : i + 3] for i in range(0, len(members) // 2, 3)]
    pairs = {(a, b): 500_000 for clique in cliques for a in clique for b in clique if a < b}
    weights = [300_000, 450_000, 500_000, 900_000]
    for _ in range(2 * len(members)):
        pairs[tuple(sorted(rng.sample(members, 2)))] = rng.choice(weights)
    if rng.random() < 0.5:
        hub, *spokes = rng.sample(members, min(len(members), 6))
        pairs.update(dict.fromkeys((tuple(sorted((hub, spoke))) for spoke in spokes), rng.choice(weights)))
    kept = [(a, b, weight) for (a, b), weight in pairs.items() if not cluster_violations([a, b], info, CONFIG)]
    clique_of = {m: number for number, clique in enumerate(cliques, first_clique) for m in clique}
    return members, kept, info, CONFIG, clique_of, TAU


def _without_names(tie):
    """What a report says whatever the ids are called (a chain is filed under its first group by name)."""
    return tie.kind, tie.strength, len(tie.candidates), tie.size if tie.kind == "tie" else None


def _reference(members, pairs, info, config, clique_of, tau):
    """The same rules written plainly and slowly, to check ``cluster_component`` against: a group is a set of ids,
    and every round sums every link again from the pairs that were not withdrawn."""
    states = {m: group_state(info[m], config.enforced_prefixes) for m in members}

    def state(*groups):
        return reduce(join_states, (states[m] for group in groups for m in group))

    def strong(weight, smaller):
        return weight >= tau * smaller if smaller == 1 else weight > tau * smaller

    group_of = {m: frozenset((m,)) for m in members}
    cliques = defaultdict(set)
    for m in members:
        if m in clique_of:
            cliques[clique_of[m]].add(m)
    for clique in map(frozenset, cliques.values()):
        if state(clique) is not None:
            group_of.update(dict.fromkeys(clique, clique))
    withdrawn_pairs, withdrawn = set(), []
    while True:
        weight = defaultdict(int)  # a link, as the set of its two groups -> the weight between them
        for a, b, pair_weight in pairs:
            if (a, b) not in withdrawn_pairs and group_of[a] != group_of[b]:
                weight[frozenset((group_of[a], group_of[b]))] += pair_weight
        strength = {
            link: Fraction(total, min(map(len, link)))
            for link, total in weight.items()
            if state(*link) is not None and strong(total, min(map(len, link)))
        }
        if not strength:
            break
        level = max(strength.values())
        strongest = [link for link in strength if strength[link] == level]
        partners = defaultdict(set)
        for one, other in strongest:
            partners[one].add(other)
            partners[other].add(one)

        def between(one, other):
            return weight.get(frozenset((one, other)), 0)

        def one_group(group, first, second):
            """Can the three be one group: the guardrails allow it, and when ``group`` has joined one of the two,
            the other still merges."""
            return state(group, first, second) is not None and (
                strong(between(group, second) + between(first, second), min(len(group) + len(first), len(second)))
                or strong(between(group, first) + between(first, second), min(len(group) + len(second), len(first)))
            )

        tied = set()
        for group, around in partners.items():
            rivals = {
                x
                for first, second in combinations(around, 2)
                if not one_group(group, first, second)
                for x in (first, second)
            }
            left = around - rivals
            if state(group, *left) is None:  # then the ones left share no family
                rivals |= {partner for partner in left if state(partner)[0] is not None}
            if rivals:
                tied |= {frozenset((group, rival)) for rival in rivals}
                withdrawn.append(("tie", level, [group, *rivals]))
        chains = []
        for link in strongest:
            if link not in tied:
                touching = [chain for chain in chains if chain & link]
                chains = [chain for chain in chains if not chain & link] + [link.union(*touching)]
        refused = [chain for chain in chains if state(*chain) is None]
        withdrawn += [("chain", level, list(chain)) for chain in refused]
        dead = tied | {link for link in strongest if any(link <= chain for chain in refused)}
        withdrawn_pairs |= {(a, b) for a, b, _weight in pairs if frozenset((group_of[a], group_of[b])) in dead}
        for chain in chains:
            if chain not in refused:
                merged = frozenset().union(*chain)
                group_of.update(dict.fromkeys(merged, merged))
    reports = set()
    for kind, level, groups in withdrawn:
        own, *others = (group_of[min(group)] for group in groups)  # where each of the groups ended
        if kind == "chain":
            own, *others = sorted({own, *others}, key=min)
        elif own in others:
            continue
        candidates = tuple(sorted({min(other) for other in others}))
        reports.add(Tie(kind, min(own), len(own), float(level) / WEIGHT_UNIT, candidates))
    return sorted(sorted(group) for group in set(group_of.values())), sorted(reports)


def test_clusters_and_reports_equal_a_plain_reference_on_random_components():
    rng = random.Random(3)
    for _ in range(1500):
        component = _random_component(rng)
        assert cluster_component(*component) == _reference(*component)


def test_clusters_ignore_id_spelling_and_input_order():
    rng = random.Random(1)
    for _ in range(300):
        members, pairs, info, _config, clique_of, _tau = component = _random_component(rng)
        clusters, ties = cluster_component(*component)
        backwards = dict(zip(sorted(members), sorted(members, reverse=True), strict=True))
        scrambled = dict(zip(members, rng.sample(members, len(members)), strict=True))
        for order in (backwards, scrambled):
            # The same evidence under other names: each id keeps its prefix, which the guardrails read, and takes
            # another id's place in the sort order. The lists come shuffled, the pairs with their ends swapped.
            new = {m: f"{m.split(':')[0]}:{order[m].replace(':', '.')}" for m in members}
            new_pairs = [(new[b], new[a], weight) for a, b, weight in rng.sample(pairs, len(pairs))]
            new_info = {new[m]: NodeInfo(curie=new[m], branches=info[m].branches, taxon=info[m].taxon) for m in members}
            new_cliques = {new[m]: clique for m, clique in clique_of.items()}
            new_members = rng.sample(sorted(new.values()), len(members))
            new_clusters, new_ties = cluster_component(new_members, new_pairs, new_info, CONFIG, new_cliques, TAU)
            assert new_clusters == sorted(sorted(new[m] for m in cluster) for cluster in clusters)
            assert sorted(map(_without_names, new_ties)) == sorted(map(_without_names, ties))


def test_an_unrelated_component_changes_nothing():
    rng = random.Random(2)
    for _ in range(300):
        one, other = _random_component(rng, "a"), _random_component(rng, "b", first_clique=100)
        clusters, ties = cluster_component(*one)
        other_clusters, other_ties = cluster_component(*other)
        (members, pairs, info, _, cliques, _), (more_members, more_pairs, more_info, _, more_cliques, _) = one, other
        both_clusters, both_ties = cluster_component(
            members + more_members, pairs + more_pairs, {**info, **more_info}, CONFIG, {**cliques, **more_cliques}, TAU
        )
        assert both_clusters == sorted(clusters + other_clusters)
        assert sorted(both_ties) == sorted(ties + other_ties)


def test_the_strongest_link_is_the_one_with_most_weight_per_id():
    one, other = ["A:1", "A:2"], ["B:1", "B:2", "MONDO:1"]
    # 1.4 between the two groups is 0.7 for each id of the smaller one: less than MONDO:2's 0.9.
    clusters, _ties = _cluster([("A:1", "B:1", 0.7), ("A:2", "B:2", 0.7), ("A:1", "MONDO:2", 0.9)], [one, other])
    assert clusters == [["A:1", "A:2", "MONDO:2"], ["B:1", "B:2", "MONDO:1"]]


def test_a_link_is_decided_again_after_its_smaller_group_grew():
    one, other = ["A:1", "A:2"], ["B:1", "B:2", "B:3"]
    # The link between the groups is 1.0 for 2 ids. Once C:1 has joined the smaller group it is 1.0 for 3 ids:
    # weaker, and still enough.
    clusters, _ties = _cluster([("A:1", "B:1", 1.0), ("A:1", "C:1", 1.0)], [one, other])
    assert clusters == [[*one, *other, "C:1"]]


def test_a_single_id_attaches_at_exactly_tau_but_two_groups_need_more():
    one, other = ["A:1", "A:2", "A:3"], ["B:1", "B:2", "B:3"]
    # Between the groups, two name matches of 0.45: 0.9 for 3 ids is exactly tau each.
    clusters, ties = _cluster([("A:1", "B:1", 0.45), ("A:2", "B:2", 0.45), ("A:3", "C:1", 0.3)], [one, other])
    assert clusters == [["A:1", "A:2", "A:3", "C:1"], ["B:1", "B:2", "B:3"]]
    assert not ties


def test_proinsulin_stays_out_of_insulin_and_of_the_proinsulin_drug():
    """The real evidence: PR:000048036 "proinsulin (human)" matches one name in insulin's Babel clique and one in the
    drug's, and those two names match each other. Had it joined insulin, the drug's 3 ids would have followed on
    0.45 + 0.45, exactly tau each, and insulin would hold a drug."""
    insulin = ["HGNC:6081", "NCBIGene:3630", "UMLS:C3813208"]
    drug = ["DRUGBANK:DB17875", "NCIT:C170355", "UMLS:C5700372"]
    name_matches = [
        ("PR:000048036", "UMLS:C3813208", 0.45),
        ("NCIT:C170355", "PR:000048036", 0.45),
        ("NCIT:C170355", "UMLS:C3813208", 0.45),
    ]
    clusters, ties = _cluster(name_matches, [insulin, drug])
    assert clusters == [sorted(drug), sorted(insulin), ["PR:000048036"]]
    assert ties == [Tie("tie", "PR:000048036", 1, 0.45, ("DRUGBANK:DB17875", "HGNC:6081"))]


def test_an_exact_tie_merges_neither_side_and_is_reported():
    # An rsid names a position, and two alleles sit at it: CAID is a one-id-per-cluster prefix. Its third link,
    # just as strong, has no rival.
    pairs = [("CAID:CA1", "DBSNP:rs1", 0.5), ("CAID:CA2", "DBSNP:rs1", 0.5), ("DBSNP:rs1", "UMLS:C1", 0.5)]
    clusters, ties = _cluster(pairs)
    assert clusters == [["CAID:CA1"], ["CAID:CA2"], ["DBSNP:rs1", "UMLS:C1"]]
    assert ties == [Tie("tie", "DBSNP:rs1", 2, 0.5, ("CAID:CA1", "CAID:CA2"))]  # the group as it ended


def test_a_tie_that_did_not_hold_is_not_reported():
    # The rsid is tied between two alleles, and a concept that is with the first one takes it there after all.
    pairs = [
        ("CAID:CA1", "DBSNP:rs1", 0.5),
        ("CAID:CA2", "DBSNP:rs1", 0.5),
        ("DBSNP:rs1", "UMLS:C1", 0.5),
        ("CAID:CA1", "UMLS:C1", 0.5),
    ]
    clusters, ties = _cluster(pairs)
    assert clusters == [["CAID:CA1", "DBSNP:rs1", "UMLS:C1"], ["CAID:CA2"]]
    assert not ties


def test_three_partners_that_only_clash_together_are_a_tie():
    # No two of the three types exclude each other, so no pair of partners is a rival: the three together are.
    categories = dict(zip(("A:1", "B:1", "C:1"), FAMILY_TRIANGLE, strict=True))
    pairs = [("A:1", "X:1", 0.5), ("B:1", "X:1", 0.5), ("C:1", "X:1", 0.5), ("D:1", "X:1", 0.5)]
    clusters, ties = _cluster(pairs, categories=categories)
    assert clusters == [["A:1"], ["B:1"], ["C:1"], ["D:1", "X:1"]]
    assert ties == [Tie("tie", "D:1", 2, 0.5, ("A:1", "B:1", "C:1"))]


def test_a_tied_link_is_withdrawn_for_good():
    # The rsid is tied between an allele's clique and another allele. Evidence that arrives with a later merge
    # (0.3 from UMLS:C1) does not add to the withdrawn 0.5: alone it is too little for two groups of 2.
    pairs = [
        ("CAID:CA1", "DBSNP:rs1", 0.5),
        ("CAID:CA2", "DBSNP:rs1", 0.5),
        ("DBSNP:rs1", "UMLS:C1", 0.45),
        ("CLINVAR:1", "UMLS:C1", 0.3),
    ]
    clusters, ties = _cluster(pairs, [["CAID:CA1", "CLINVAR:1"]])
    assert clusters == [["CAID:CA1", "CLINVAR:1"], ["CAID:CA2"], ["DBSNP:rs1", "UMLS:C1"]]
    assert len(ties) == 1


def test_equal_links_to_groups_that_can_be_one_group_are_all_applied():
    one, other = ["A:1", "A:2", "A:3"], ["B:1", "B:2", "B:3"]
    # X:1's 0.45 is too little for 3 ids, but the 0.5 between the two groups makes up for it: 0.95 for 3 ids.
    clusters, ties = _cluster([("A:1", "X:1", 0.45), ("B:1", "X:1", 0.45), ("A:2", "B:2", 0.5)], [one, other])
    assert clusters == [[*one, *other, "X:1"]] and not ties
    # Nothing has to be made up for a single id, or for groups the links are strong enough for: 0.9 for 3 ids each.
    clusters, ties = _cluster([("A:1", "X:1", 0.45), ("C:1", "X:1", 0.45)], [one])
    assert clusters == [[*one, "C:1", "X:1"]] and not ties
    clusters, ties = _cluster([("A:1", "X:1", 1.0), ("B:1", "X:1", 1.0)], [one, other])
    assert clusters == [[*one, *other, "X:1"]] and not ties


def test_a_valid_babel_clique_is_never_split():
    # Each half of the clique has stronger evidence for another disease than the clique's own 0.5.
    clusters, _ties = _cluster([("MESH:D1", "MONDO:1", 1.0), ("MONDO:2", "UMLS:C1", 1.0)], [["MESH:D1", "UMLS:C1"]])
    assert ["MESH:D1", "UMLS:C1"] in clusters


def test_a_babel_clique_that_breaks_a_guardrail_does_not_start_as_one_group():
    taxa = {"NCBIGene:1": HUMAN, "NCBIGene:2": MOUSE}
    clusters, _ties = _cluster([], [["NCBIGene:1", "NCBIGene:2", "UMLS:C1"]], taxa)
    assert clusters == [["NCBIGene:1"], ["NCBIGene:2"], ["UMLS:C1"]]


def test_a_guardrail_blocks_a_merge_the_pairwise_prune_cannot_see():
    taxa = {"NCBIGene:1": HUMAN, "NCBIGene:2": MOUSE}
    pairs = [("NCBIGene:1", "UMLS:C1", 1.0), ("NCBIGene:2", "UMLS:C1", 0.5), ("UMLS:C1", "UMLS:C2", 0.5)]
    clusters, ties = _cluster(pairs, taxa=taxa)
    assert clusters == [["NCBIGene:1", "UMLS:C1", "UMLS:C2"], ["NCBIGene:2"]]
    assert not ties  # the forbidden link is no rival to the equally strong one that is allowed


def test_a_chain_of_equal_links_the_guardrails_forbid_merges_nothing():
    taxa = {"NCBIGene:1": HUMAN, "NCBIGene:2": MOUSE}
    chain = [("NCBIGene:1", "UMLS:C1", 0.5), ("UMLS:C1", "UMLS:C2", 0.5), ("NCBIGene:2", "UMLS:C2", 0.5)]
    # ...for good: when UMLS:C1 later joins another group on a weaker link, the chain's links are not tried again.
    clusters, ties = _cluster([*chain, ("UMLS:C1", "UMLS:X1", 0.4)], [["UMLS:X1", "UMLS:X2"]], taxa)
    assert clusters == [["NCBIGene:1"], ["NCBIGene:2"], ["UMLS:C1", "UMLS:X1", "UMLS:X2"], ["UMLS:C2"]]
    # No group had to choose, so it is no tie: the chain is reported as one, its groups named as they ended.
    assert ties == [Tie("chain", "NCBIGene:1", 1, 0.5, ("NCBIGene:2", "UMLS:C1", "UMLS:C2"))]


def test_two_refused_chains_that_end_as_the_same_groups_are_one_row():
    taxa = {"NCBIGene:1": HUMAN, "NCBIGene:2": MOUSE, "NCBIGene:3": HUMAN, "NCBIGene:4": MOUSE}
    chains = [(a, b, 0.5) for a, b in [("NCBIGene:1", "UMLS:C1"), ("UMLS:C1", "UMLS:C2"), ("NCBIGene:2", "UMLS:C2")]]
    chains += [(a, b, 0.5) for a, b in [("NCBIGene:3", "UMLS:C3"), ("UMLS:C3", "UMLS:C4"), ("NCBIGene:4", "UMLS:C4")]]
    # Weaker links then join the two chains' groups in pairs, so the two reports name the same four groups.
    joins = [("NCBIGene:1", "NCBIGene:3"), ("UMLS:C1", "UMLS:C3"), ("UMLS:C2", "UMLS:C4"), ("NCBIGene:2", "NCBIGene:4")]
    _clusters, ties = _cluster([*chains, *((a, b, 0.4) for a, b in joins)], taxa=taxa)
    assert ties == [Tie("chain", "NCBIGene:1", 2, 0.5, ("NCBIGene:2", "UMLS:C1", "UMLS:C2"))]


def test_the_tie_report_reads_back_as_it_was_written(tmp_path):
    # Strengths one millionth apart, and an id with a space in it.
    ties = [
        Tie("tie", "DBSNP:rs1", 2, 1.000001, ("CAID:CA1", "CAID:CA2")),
        Tie("tie", "DBSNP:rs1", 2, 1.000002, ("CAID:CA1", "HCPCS:Level 1: A0000-A9999")),
    ]
    log_ties(ties, tmp_path / "er_ties.tsv")
    _header, *rows = (tmp_path / "er_ties.tsv").read_text().splitlines()
    fields = [row.split("\t") for row in rows]
    assert [
        (kind, group, int(size), float(strength), tuple(json.loads(c))) for kind, group, size, strength, c in fields
    ] == ties
