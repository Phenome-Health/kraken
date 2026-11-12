"""Last-resort category and taxon guesses from an identifier's prefix.

Entity resolution types and taxons every identifier from, in order: Babel's own facts for it (see id_facts), the
categories/taxon a source gave it, and only then these prefix guesses. A prefix heuristic must never override a real
answer, so callers apply these AFTER source values -- never inside a facts lookup.
"""

from kraken.utils.constants import ORGANISM_TAXON_CATEGORY

# Unambiguous prefix -> Biolink category. Conservative on purpose; extend only as safe.
PREFIX_CATEGORY: dict[str, str] = {
    "HGNC": "biolink:Gene",
    "NCBIGene": "biolink:Gene",
    "UniProtKB": "biolink:Protein",
    "PR": "biolink:Protein",
    "NCBITaxon": ORGANISM_TAXON_CATEGORY,
    "CHEBI": "biolink:ChemicalEntity",
    "PUBCHEM.COMPOUND": "biolink:SmallMolecule",
    "KEGG.COMPOUND": "biolink:SmallMolecule",
    "HMDB": "biolink:SmallMolecule",
    "INCHIKEY": "biolink:ChemicalEntity",
    "MONDO": "biolink:Disease",
    # HP terms are usually phenotypes but sometimes diseases, so back off to the shared parent.
    "HP": "biolink:DiseaseOrPhenotypicFeature",
    "UBERON": "biolink:AnatomicalEntity",
    "CL": "biolink:Cell",
    # NOTE: DRUGBANK is intentionally NOT here -- it spans small molecules AND biologics/protein drugs, so no
    # single category is safe.
    # NOTE: REACT, ENSEMBL, bare GO and KEGG are intentionally NOT here either: each spans several kinds of thing
    # (Reactome R-HSA ids are pathways, reactions AND physical entities like modified-protein states), so a
    # blanket guess would mis-type many nodes.
}

# Prefix -> taxon for single-species nomenclature authorities (the prefix itself DEFINES the species), at species
# rank like every other taxon in the graph. Multi-species prefixes (NCBIGene, UniProtKB, ENSEMBL, Xenbase) are
# deliberately omitted.
TAXON_BY_PREFIX: dict[str, str] = {
    "HGNC": "NCBITaxon:9606",  # human
    "MGI": "NCBITaxon:10090",  # mouse
    "RGD": "NCBITaxon:10116",  # rat -- a default only, see DEFAULT_ONLY_TAXON_PREFIXES
    "ZFIN": "NCBITaxon:7955",  # zebrafish
    "FB": "NCBITaxon:7227",  # fruit fly (FlyBase)
    "FlyBase": "NCBITaxon:7227",
    "WB": "NCBITaxon:6239",  # C. elegans (WormBase)
    "WormBase": "NCBITaxon:6239",
    "SGD": "NCBITaxon:4932",  # S. cerevisiae (the species, not SGD's reference strain S288C)
    "PomBase": "NCBITaxon:4896",  # S. pombe
    "dictyBase": "NCBITaxon:44689",  # D. discoideum
    "TAIR": "NCBITaxon:3702",  # A. thaliana
}

# Prefixes above that do NOT define the species. RGD names dog and pig genes as well as rat ones, so rat is only its
# default: entity resolution gives an RGD id the species of its Babel clique where that is known to be another one
# (build._clique_taxon_overrides). RGD stays in the table for the rat ids nothing else gives a species. Every other
# prefix's guess agreed with all of its Babel cliques in the 2.4.0-dev3 build.
DEFAULT_ONLY_TAXON_PREFIXES: frozenset[str] = frozenset({"RGD"})


def infer_category(curie: str) -> str | None:
    return PREFIX_CATEGORY.get(curie.split(":", 1)[0])


def infer_taxon(curie: str) -> str | None:
    """Backup taxon from a single-species nomenclature prefix (or None)."""
    return TAXON_BY_PREFIX.get(curie.split(":", 1)[0])
