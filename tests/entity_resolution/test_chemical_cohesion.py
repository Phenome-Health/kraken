"""Putting back the chemical Babel cliques that clustering split."""

import pytest

from kraken.entity_resolution.chemical_cohesion import plan_chemical_clique_repairs
from kraken.entity_resolution.families import ALL_FAMILIES
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo

CHEMICAL = frozenset({"chemical"})
GENE_PROTEIN = frozenset({"gene_protein"})
SMALL_MOLECULE = ("biolink:SmallMolecule",)
PROTEIN = ("biolink:Protein",)
GENE = ("biolink:Gene",)
FAMILY = {"biolink:SmallMolecule": CHEMICAL, "biolink:Protein": GENE_PROTEIN, "biolink:Gene": GENE_PROTEIN}


def _plan(tmp_path, clusters, cliques, categories=None, *, names=None, pairs=None, outliers=()):
    """Every id is a SmallMolecule unless ``categories`` says otherwise. Without ``pairs``, every two ids of a cluster
    are a pair."""
    categories = categories or {}
    path = tmp_path / "cliques.tsv"
    path.write_text("".join(f"{m}\t{c[0]}\t{len(c)}\n" for c in cliques for m in c))
    if pairs is None:
        by_cluster = {}
        for curie, cluster in sorted(clusters.items()):
            by_cluster.setdefault(cluster, []).append(curie)
        pairs = [(a, b) for ids in by_cluster.values() for a in ids for b in ids if a < b]

    def info_of(curie):
        cats = categories.get(curie, SMALL_MOLECULE)
        return NodeInfo(curie=curie, branches=FAMILY[cats[0]] if cats else ALL_FAMILIES)

    typed = {curie: categories.get(curie, SMALL_MOLECULE) for curie in clusters}
    return plan_chemical_clique_repairs(
        clusters,
        path,
        typed,
        info_of,
        (names or {}).get,
        lambda ids: [(a, b) for a, b in pairs if a in ids and b in ids],
        GuardrailConfig(),
        100,
        outliers,
    )


def _sorted(clusters):
    return sorted(sorted(c) for c in clusters)


def _unchanged(plan):
    return not plan.replaced and not plan.remade and not plan.retyped


NAC = "N-acetyl-S-allylcysteine"


def test_a_protein_typed_concept_joins_the_chemical_it_names(tmp_path):
    """N-acetyl-S-allylcysteine: UMLS T116 puts it in Babel's Protein compendium, with its MeSH id."""
    clusters = {"UMLS:C0526228": 0, "MESH:C054604": 1, "CHEBI:189718": 1, "INCHIKEY:LKRAEHUDIUJBSF-ZETCQYMHSA-N": 1}
    cliques = [
        ["UMLS:C0526228", "MESH:C054604"],
        ["CHEBI:189718", "MESH:C054604", "INCHIKEY:LKRAEHUDIUJBSF-ZETCQYMHSA-N"],
    ]
    names = dict.fromkeys(["UMLS:C0526228", "MESH:C054604", "CHEBI:189718"], NAC)
    plan = _plan(tmp_path, clusters, cliques, {"UMLS:C0526228": PROTEIN}, names=names)
    assert plan.retyped == {"UMLS:C0526228": SMALL_MOLECULE}
    assert _sorted(plan.remade) == [sorted(clusters)]


@pytest.mark.parametrize(
    ("proteins", "link", "chemical", "names"),
    [
        # An antibody-drug conjugate and MeSH's "SC002", whose SmallMolecule clique is a cephalosporin.
        (
            ["NCIT:C123284", "NCIT:C203531", "UMLS:C0967624"],
            "MESH:C430544",
            ["PUBCHEM.COMPOUND:136693871", "PUBCHEM.COMPOUND:177582760"],
            {
                "NCIT:C123284": "Anti-DLL3 Antibody-drug Conjugate SC-002",
                "NCIT:C203531": "Anti-DLL3 ADC ZL-1310",
                "UMLS:C0967624": "SC002",
                "PUBCHEM.COMPOUND:136693871": "SC-002",
            },
        ),
        # The family of growth factors, and MeSH's plural, whose SmallMolecule clique is one peptide's structure.
        (
            ["NCIT:C20426", "UMLS:C0016026"],
            "MESH:D005346",
            ["PUBCHEM.COMPOUND:133126676", "INCHIKEY:MSJSNXQCIRMHFJ-WIRNGOBFSA-N"],
            {
                "NCIT:C20426": "Fibroblast Growth Factor",
                "UMLS:C0016026": "Fibroblast Growth Factor",
                "PUBCHEM.COMPOUND:133126676": "Fibroblast Growth Factors",
            },
        ),
    ],
)
def test_a_concept_that_names_no_member_of_the_chemical_stays_a_protein(tmp_path, proteins, link, chemical, names):
    clusters = dict.fromkeys(proteins, 0) | dict.fromkeys([link, *chemical], 1)
    plan = _plan(
        tmp_path, clusters, [[*proteins, link], [link, *chemical]], dict.fromkeys(proteins, PROTEIN), names=names
    )
    assert _unchanged(plan)


@pytest.mark.parametrize(
    ("protein", "link", "pubchem", "name", "others"),
    [
        # PubChem titles ergothioneine "Thioredoxins", and a chlorinated methanoindene "Peptide Elongation Factor 2":
        # the MeSH names linked to them.
        (
            "UMLS:C0039938",
            "MESH:D013879",
            "PUBCHEM.COMPOUND:3037043",
            "Thioredoxins",
            {"CHEBI:134344": "ergothioneine(1+)", "INCHIKEY:SSISHJJTAXXQAX-ZETCQYMHSA-O": None},
        ),
        (
            "UMLS:C0059037",
            "MESH:D020652",
            "PUBCHEM.COMPOUND:148754",
            "Peptide Elongation Factor 2",
            {"INCHIKEY:RYXGAWJESREXQK-UHFFFAOYSA-N": None},
        ),
    ],
)
def test_a_name_only_pubchem_gives_the_chemical_moves_nothing(tmp_path, protein, link, pubchem, name, others):
    clusters = {protein: 0} | dict.fromkeys([link, pubchem, *others], 1)
    names = dict.fromkeys([protein, link, pubchem], name) | others
    cliques = [[protein, link], [link, pubchem, *others]]
    assert _unchanged(_plan(tmp_path, clusters, cliques, {protein: PROTEIN}, names=names))


@pytest.mark.parametrize(
    ("protein", "links", "chemical", "name"),
    [
        # Methotrexate's brands sit with "Methopterin": DrugBank's methopterin is the link.
        ("UMLS:C5441719", ["DRUGBANK:DB16620"], ["PUBCHEM.COMPOUND:135434845"], "Methopterin"),
        # Tocilizumab: the node those links reach also holds DRUGBANK:DB08208, a small molecule.
        ("MESH:C502936", ["DRUGBANK:DB06273", "UMLS:C1609165"], ["CHEBI:64360", "DRUGBANK:DB08208"], "tocilizumab"),
        # Asparaginase, the enzyme: its links sit in a node with asparaginate.
        ("MESH:D001215", ["DRUGBANK:DB00023", "UMLS:C0003993"], ["CHEMBL.COMPOUND:CHEMBL2108989"], "asparaginase"),
    ],
)
def test_a_link_that_is_not_a_mesh_id_moves_nothing(tmp_path, protein, links, chemical, name):
    clusters = {protein: 0} | dict.fromkeys(links + chemical, 1)
    cliques = [[protein, *links], [*links, *chemical]]
    names = {protein: name, chemical[0]: name}
    assert _unchanged(_plan(tmp_path, clusters, cliques, {protein: PROTEIN}, names=names))


def test_a_concept_in_a_cluster_with_a_protein_database_id_stays_a_protein(tmp_path):
    """UMLS's "Erythropoietin" shares MESH:D004921 with a PubChem record, but it belongs with PR:000007141."""
    clusters = {"UMLS:C0014822": 0, "PR:000007141": 0, "MESH:D004921": 1, "PUBCHEM.COMPOUND:92043599": 1}
    cliques = [["UMLS:C0014822", "MESH:D004921"], ["PUBCHEM.COMPOUND:92043599", "MESH:D004921"]]
    categories = {"UMLS:C0014822": PROTEIN, "PR:000007141": PROTEIN}
    names = dict.fromkeys(["UMLS:C0014822", "PUBCHEM.COMPOUND:92043599"], "Erythropoietin")
    assert _unchanged(_plan(tmp_path, clusters, cliques, categories, names=names))


def test_a_member_also_in_a_gene_clique_stays_a_protein(tmp_path):
    clusters = {"UMLS:C1": 0, "MESH:D1": 1, "UNII:1": 1, "NCBIGene:1": 2}
    cliques = [["UMLS:C1", "MESH:D1"], ["UNII:1", "MESH:D1"], ["NCBIGene:1", "UMLS:C1"]]
    categories = {"UMLS:C1": PROTEIN, "NCBIGene:1": GENE}
    assert _unchanged(_plan(tmp_path, clusters, cliques, categories, names={"UMLS:C1": "x1", "UNII:1": "x1"}))


@pytest.mark.parametrize(
    ("proteins", "link", "chemical", "name", "beside", "pairs"),
    [
        # S-adenosylhomocysteine: LOINC's part is there by its equivalence to UMLS's concept, and would be left alone.
        (
            ["NCIT:C103149", "UMLS:C0036001"],
            "MESH:D012435",
            "CHEBI:16680",
            "S-Adenosylhomocysteine",
            ["LOINC:LP200579-3"],
            [("LOINC:LP200579-3", "UMLS:C0036001")],
        ),
        # Ristocetin: LOINC's part likewise, with SNOMED's concept tied to it, so what stays holds together.
        (
            ["NCIT:C152234", "UMLS:C0035654"],
            "MESH:D012310",
            "CHEBI:85129",
            "Ristocetin",
            ["LOINC:LP34424-9", "SNOMEDCT:31827005"],
            [("LOINC:LP34424-9", "UMLS:C0035654"), ("LOINC:LP34424-9", "SNOMEDCT:31827005")],
        ),
    ],
)
def test_a_concept_with_ids_of_its_own_beside_it_stays(tmp_path, proteins, link, chemical, name, beside, pairs):
    clusters = dict.fromkeys([*proteins, *beside], 0) | dict.fromkeys([link, chemical], 1)
    pairs = [*pairs, (proteins[0], proteins[1]), (chemical, link)]
    categories = dict.fromkeys(proteins, PROTEIN) | dict.fromkeys(beside, ())
    names = dict.fromkeys([*proteins, chemical], name)
    plan = _plan(tmp_path, clusters, [[*proteins, link], [link, chemical]], categories, names=names, pairs=pairs)
    assert _unchanged(plan)


# Made-up InChIKeys keep the real shape: the block before "-" is the skeleton, shared by stereo forms.
CORTISONE = ["CHEBI:16962", "INCHIKEY:C-Z", "HMDB:HMDB0002802", "CHEBI:14026", "MESH:D003348", "UMLS:C0010137"]
NO_STEREO = ["PUBCHEM.COMPOUND:521440", "INCHIKEY:C-U"]


def test_most_of_a_clique_with_another_inchikey_joins_its_own_and_nothing_else_moves(tmp_path):
    """Cortisone: most of its clique sits with the no-stereo InChIKey; RefMet's id stays with the InChIKey it lists,
    and LOINC's "Cortisone", outside the clique, stays where label propagation put it."""
    clusters = dict.fromkeys([*CORTISONE[:2], "RM:0128363"], 0) | dict.fromkeys([*CORTISONE[2:], *NO_STEREO], 1)
    clusters["LOINC:LP18544-4"] = 1
    pairs = [(a, b) for a in CORTISONE for b in CORTISONE if a < b] + [("INCHIKEY:C-Z", "RM:0128363")]
    pairs += [("PUBCHEM.COMPOUND:521440", "INCHIKEY:C-U"), ("LOINC:LP18544-4", "PUBCHEM.COMPOUND:521440")]
    pairs += [("LOINC:LP18544-4", "UMLS:C0010137")]
    plan = _plan(tmp_path, clusters, [CORTISONE, NO_STEREO], pairs=pairs)
    assert _sorted(plan.remade) == _sorted([[*CORTISONE, "RM:0128363"], [*NO_STEREO, "LOINC:LP18544-4"]])


def test_most_of_a_clique_with_another_skeleton_is_left_alone(tmp_path):
    """MeSH's GM3 ganglioside, a class, sits with one GM3 species and is in the clique of another, GM3(d18:1/24:1(15Z)):
    two compounds, not two stereo forms of one."""
    clique = ["INCHIKEY:PFJKOHUKELZMLE-A", "MESH:D005679", "UMLS:C0016897"]
    clusters = {clique[0]: 0, "LM:SP0601AJ07": 0} | dict.fromkeys(
        [*clique[1:], "INCHIKEY:NVMUTWNBCCAIGX-A", "LM:GM3"], 1
    )
    assert _unchanged(_plan(tmp_path, clusters, [clique]))


IBUPROFEN_SMILES = "SMILES:CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1"


@pytest.mark.parametrize(
    ("clique", "stays", "ties"),
    [
        # Ibuprofen: the SMILES in the donor is the structure of the clique's InChIKey, tied to the members that would
        # leave and to nothing that stays.
        (
            ["INCHIKEY:I-R", "CHEBI:47835", "CHEMBL.COMPOUND:CHEMBL427526", "UNII:99W8H60N62"],
            [IBUPROFEN_SMILES, "PUBCHEM.COMPOUND:3672", "INCHIKEY:I-U"],
            [
                (IBUPROFEN_SMILES, "CHEBI:47835"),
                ("PUBCHEM.COMPOUND:3672", "INCHIKEY:I-U"),
                ("CHEBI:47835", "PUBCHEM.COMPOUND:3672"),
            ],
        ),
        # Mycothione: KEGG's record is what ties another InChIKey to the clique, and would leave it alone.
        (
            ["INCHIKEY:M-M", "CHEBI:16086", "CHEBI:14627", "CHEBI:7041", "KEGG.COMPOUND:C11502"],
            ["INCHIKEY:M-K"],
            [("INCHIKEY:M-K", "KEGG.COMPOUND:C11502")],
        ),
    ],
)
def test_a_move_that_leaves_the_donor_in_pieces_is_not_made(tmp_path, clique, stays, ties):
    clusters = {clique[0]: 0} | dict.fromkeys([*clique[1:], *stays], 1)
    pairs = [(a, b) for a in clique for b in clique if a < b] + ties
    assert _unchanged(_plan(tmp_path, clusters, [clique], pairs=pairs))


def test_a_move_that_leaves_a_smiles_behind_is_not_made(tmp_path):
    """All-trans-neoxanthin: LIPID MAPS lists ChEBI's InChIKey with a SMILES whose InChIKey is that of the PubChem
    clique, which would leave the SMILES with the other stereoisomer."""
    clique = ["INCHIKEY:PGYAYSRVSAJXTE-CLONMANBSA-N", "PUBCHEM.COMPOUND:5281247", "CAS:30743-41-0"]
    smiles = (
        "SMILES:C/C(C=C=C1C(C)(C)C[C@H](O)C[C@@]1(C)O)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(C)\\C=C\\"
        "[C@@]12O[C@]1(C)C[C@@H](O)CC2(C)C"
    )
    chebi = ["CHEBI:32446", "INCHIKEY:PGYAYSRVSAJXTE-MTYISEJWSA-N"]
    clusters = {clique[0]: 0} | dict.fromkeys([*clique[1:], *chebi, "LM:PR01070290", smiles], 1)
    assert _unchanged(_plan(tmp_path, clusters, [clique, chebi]))


@pytest.mark.parametrize(
    "away_with",
    [
        # Alloaromadendrene: ChEBI's "(+)-alloaromadendrene" is in the clique of the (-) form's InChIKey.
        [],
        # Two RefMet ids, one for each piece: apart because of RefMet's duplicate, not two structures.
        ["RM:0000002"],
    ],
)
def test_a_clique_with_no_other_inchikey_in_its_clusters_is_left_alone(tmp_path, away_with):
    clique = ["INCHIKEY:ITYNGVSTWVVPIC-DHGKCCLASA-N", "PUBCHEM.COMPOUND:10899740", "CHEBI:228351", "MESH:C419726"]
    clique += ["UMLS:C0963965", "CHEMBL.COMPOUND:CHEMBL3120657"]
    clusters = dict.fromkeys([*clique[:2], "LM:PR0103430003", "RM:0000001"], 0) | dict.fromkeys(clique[2:], 1)
    clusters |= dict.fromkeys(away_with, 1)
    assert _unchanged(_plan(tmp_path, clusters, [clique]))


def test_a_few_ids_that_left_their_clique_for_another_structure_stay_put(tmp_path):
    """GtoPdb's "L-tryptophan" sits in Babel's no-stereo clique; its name is what puts it with L-tryptophan."""
    l_form = ["CHEBI:16828", "INCHIKEY:W-L"]
    no_stereo = ["CHEBI:27897", "INCHIKEY:W-U", "HMDB:HMDB0030396", "GTOPDB:717"]
    clusters = dict.fromkeys([*l_form, "GTOPDB:717"], 0) | dict.fromkeys(no_stereo[:3], 1)
    assert _unchanged(_plan(tmp_path, clusters, [l_form, no_stereo]))


def test_a_clique_split_evenly_with_other_structures_is_left_alone(tmp_path):
    """Xanthatin: GtoPdb's "xanthatin" and ChEBI's "8-epixanthatin" are half of Babel's no-stereo clique, each with its
    own stereoisomer by name."""
    no_stereo = ["INCHIKEY:X-U", "PUBCHEM.COMPOUND:282783", "GTOPDB:13700", "CHEBI:180746"]
    xanthatin = ["CHEBI:10058", "INCHIKEY:X-X"]
    epixanthatin = ["HMDB:HMDB0036473", "INCHIKEY:X-E"]
    clusters = dict.fromkeys(no_stereo[:2], 0) | {"GTOPDB:13700": 1, "CHEBI:180746": 2}
    clusters |= dict.fromkeys(xanthatin, 1) | dict.fromkeys(epixanthatin, 2)
    assert _unchanged(_plan(tmp_path, clusters, [no_stereo, xanthatin, epixanthatin]))


def test_a_babel_outlier_is_not_a_member_of_its_clique(tmp_path):
    """babel_outliers took tropacocaine out of benzoyltropeine's clique. Counted in, it and the one other member away
    would be most of the clique, sitting with another InChIKey, and would be moved."""
    clique = ["INCHIKEY:B-B", "CHEBI:3043", "CHEBI:9755"]
    tropacocaine = ["INCHIKEY:B-T", "UNII:1I92X32F6H"]
    clusters = {"INCHIKEY:B-B": 0} | dict.fromkeys([*clique[1:], *tropacocaine], 1)
    pairs = [("CHEBI:3043", "CHEBI:9755"), ("CHEBI:9755", "UNII:1I92X32F6H"), tuple(tropacocaine)]
    plan = _plan(tmp_path, clusters, [clique, tropacocaine], pairs=pairs, outliers={"CHEBI:9755"})
    assert _unchanged(plan)
