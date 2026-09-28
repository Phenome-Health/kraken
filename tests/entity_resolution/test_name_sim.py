"""Tests for name normalization and name-similarity grouping."""

from kraken.entity_resolution.families import ALL_FAMILIES
from kraken.entity_resolution.name_sim import (
    group_by_normalized_name,
    is_droppable,
    name_keys,
    name_similarity_edges,
    normalize_name,
)

CHEMICAL = frozenset({"chemical"})
DISEASE = frozenset({"disease_pheno"})


def same_key(a: str, b: str, branches: frozenset[str] = ALL_FAMILIES, other: frozenset[str] | None = None) -> bool:
    """Whether the two names match: ``a`` on a node of ``branches``, ``b`` on one of ``other`` (default: the same)."""
    return bool(name_keys(a, branches) & name_keys(b, branches if other is None else other))


def test_name_key_ignores_spacing_hyphens_and_possessives():
    assert same_key("LY-2940094", "LY2940094")
    assert same_key("AM 404", "AM404")
    assert same_key("12(R)-HETE", "12R-HETE")
    assert same_key("Vitamin B 12", "Vitamin B12")
    assert same_key("Parkinson's disease", "Parkinson disease")
    assert same_key("Alzheimer’s disease 10", "Alzheimer disease 10")


def test_name_key_keeps_what_distinguishes_the_entity():
    # stereo signs: the enantiomers, and either vs the racemate or an unsigned name
    assert not same_key("EPICHLOROHYDRIN, (+)-", "EPICHLOROHYDRIN, (-)-")
    assert not same_key("(+)-camphor", "camphor")
    assert same_key("VORICONAZOLE, (+/-)-", "voriconazole (±)")
    # charges and primes
    assert not same_key("FMNH(.)(2-)", "FMNH2")
    assert not same_key("Myricetin 3'-glucoside", "Myricetin 3-glucoside")
    assert not same_key("3',6-Disinapoylsucrose", "3,6'-Disinapoyl sucrose")
    assert same_key("2'-O-Methyl Uridine", "2′-O-methyluridine")
    # the separator between two numbers: that there is one, a decimal point, and lipid sn-positions known ("/") or not
    assert not same_key("1,2-diacylglycerol", "12-diacylglycerol")
    assert not same_key("0.1 [answer]", "0-1 [answer]")
    assert not same_key("PC 14:0/20:0", "PC 14:0_20:0")
    assert same_key("TG(17:0/17:0/19:0)", "TG 17:0/17:0/19:0")
    # ...but a locant break is one however it is written
    assert same_key("14,15-EpETE", "14(15)-EpETE")
    assert same_key("1,1,1-trichloroethane", "1-1-1-trichloroethane")
    assert same_key("A + O2 => (3S)-B", "A + O2 <=> (3S)-B")


def test_name_key_singularizes_outside_chemistry_only():
    assert same_key("Jejunal Neoplasms", "jejunal neoplasm", DISEASE)
    assert same_key("Allergies", "allergy", DISEASE)
    assert same_key("status epilepticus", "Status Epilepticus", DISEASE)  # not a plural
    # chemistry names a CLASS with the plural; a wildcard node may be a chemical
    assert not same_key("uridines", "uridine", CHEMICAL)
    assert not same_key("Resorcinols", "resorcinol", ALL_FAMILIES)
    # a typed plural meets the singular on any node, and identical names always match
    assert same_key("Retinal Diseases", "retinal disease", DISEASE, ALL_FAMILIES)
    finding = "Physical findings.mandibular condyle"
    assert same_key(finding, finding, ALL_FAMILIES, DISEASE)


def test_name_key_keeps_spacing_for_organisms_genes_and_variants():
    # distinct NCBI taxa differ only by punctuation
    organism = frozenset({"organism"})
    assert not same_key("Burkholderia sp. S2", "Burkholderia sp. S-2", organism)
    assert not same_key("Actias", "Actia", organism)
    assert name_keys("Homo sapiens", organism) == {normalize_name("Homo sapiens")}
    assert name_keys("IL-6", frozenset({"gene_protein"})) == {"il 6"}
    # an untyped node still meets them
    assert same_key("Mobiluncus sp", "Mobiluncus sp.", ALL_FAMILIES, organism)


def test_normalize_basic():
    assert normalize_name("  Adams-Oliver Syndrome 1 ") == "adams oliver syndrome 1"
    assert normalize_name("Café-Résumé") == "cafe resume"  # combining accents stripped
    assert normalize_name("β-amyloid") == "β amyloid"  # non-combining Greek letter kept as-is
    assert normalize_name(None) == ""
    assert normalize_name("HbA1c!!") == "hba1c"


def test_is_droppable():
    assert is_droppable("")
    assert is_droppable("ab", min_length=3)
    assert is_droppable("123")  # purely numeric
    assert is_droppable("point in time")  # stoplist
    assert not is_droppable("insulin")


def test_group_primary_names_only():
    pairs = [
        ("A:1", "Insulin"),
        ("B:1", "insulin"),
        ("C:1", "INSULIN"),
        ("D:1", "glucose"),
        ("E:1", "12345"),  # numeric -> dropped
        ("F:1", "ab"),  # too short -> dropped
    ]
    groups = group_by_normalized_name(pairs)
    assert groups == {"insulin": ["A:1", "B:1", "C:1"]}  # glucose singleton dropped


def test_name_similarity_edges_and_cap():
    groups = {"insulin": ["A:1", "B:1", "C:1"]}
    edges = list(name_similarity_edges(groups))
    assert set(edges) == {("A:1", "B:1"), ("A:1", "C:1"), ("B:1", "C:1")}
    # oversized group skipped
    big = {"x": [f"N:{i}" for i in range(50)]}
    assert list(name_similarity_edges(big, group_cap=40)) == []


def test_an_rsid_is_an_identifier_not_a_name():
    """5.1M CAID and 5.0M DBSNP nodes are named "rs10154897" and the like -- an rsid names the POSITION, so every
    allele there shares it. Matching on that is identifier matching: it glued every allele at a position together,
    and the CAID guardrail then split them all apart again (which is what filled the build log). The real
    allele/position link is carried as `member_of` edges."""
    from kraken.entity_resolution.name_sim import is_droppable, normalize_name

    assert is_droppable(normalize_name("rs10154897"))
    assert is_droppable(normalize_name("RS10154897"))
    assert not is_droppable(normalize_name("rs1 variant of BRCA1"))  # a real name that merely mentions one
    assert not is_droppable(normalize_name("metformin"))
