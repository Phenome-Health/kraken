"""
Species whose genes the NCBI Gene harmonizer ingests.

NCBI's all-species `gene_info` holds ~72M genes across ~54k taxa -- several times the size of the rest of
KRAKEN combined, and overwhelmingly from genomes of no biomedical interest. This list scopes the ingest to
organisms that matter for a biomedical/wellness knowledge graph.

WHY THIS IS CURATED RATHER THAN DERIVED. Two automatic approaches were measured and both fail:

  - Ranking taxa by their gene count in gene_info gives 65% coverage of the NCBIGene ids already in KRAKEN,
    for 8.0M genes. Its top entries are a corroboree frog, wheat, an axolotl and a banana -- it ranks genome
    size and sequencing effort, not relevance. Rat doesn't even make the top 100.
  - Ranking by how many genes carry a dbXref does better (88%, 5.4M genes) but still leads with Trichomonas,
    locusts and a cattle tick.

`gene_info` simply carries no signal about biomedical relevance, so no ranking over it can supply one. The
remaining automatic option -- deriving the list from the taxa KRAKEN already contains -- is circular: it can
never introduce an organism the graph doesn't yet have, freezing coverage at today's sources. Curation is
what lets the list be both independent and forward-looking.

For reference, this list covers 98.4% of the NCBIGene ids in KRAKEN 2.1.0 while ingesting ~2.3M genes.
Re-measure with scripts/audit_ncbigene_taxon_allowlist.py, which reports coverage against a built node set
and names the biggest uncovered species -- use it to decide what to ADD here; it never rewrites this file.

RANK. Entries may be written at whatever rank names the organism (dog is published as the subspecies
Canis lupus familiaris, 9615). The harmonizer rolls both these entries and each gene's taxon up to species
before comparing, so no entry needs to be species-rank itself.
"""

# Curated by organism group. Each entry is an NCBI tax_id; the comment is its scientific name.
TAXON_ALLOWLIST: frozenset[str] = frozenset(
    {
        # --- Mammals ---
        "9606",  # Homo sapiens
        "10090",  # Mus musculus
        "10116",  # Rattus norvegicus
        "9913",  # Bos taurus
        "9823",  # Sus scrofa
        "9615",  # Canis lupus familiaris
        "9685",  # Felis catus
        "9796",  # Equus caballus
        "9940",  # Ovis aries
        "9925",  # Capra hircus
        "9986",  # Oryctolagus cuniculus
        "9544",  # Macaca mulatta
        "9598",  # Pan troglodytes
        "9541",  # Macaca fascicularis
        "10141",  # Cavia porcellus
        "10029",  # Cricetulus griseus (CHO -- the workhorse cell line for biologics manufacturing)
        "9361",  # Dasypus novemcinctus
        "13616",  # Monodelphis domestica
        # --- Other vertebrates ---
        "9031",  # Gallus gallus
        "7955",  # Danio rerio
        "8364",  # Xenopus tropicalis
        "8355",  # Xenopus laevis
        "8090",  # Oryzias latipes
        "69293",  # Gasterosteus aculeatus
        "9103",  # Meleagris gallopavo
        "59729",  # Taeniopygia guttata
        "28377",  # Anolis carolinensis
        # --- Invertebrates ---
        "7227",  # Drosophila melanogaster
        "6239",  # Caenorhabditis elegans
        "7460",  # Apis mellifera
        "7165",  # Anopheles gambiae
        "7159",  # Aedes aegypti
        "7668",  # Strongylocentrotus purpuratus
        "6183",  # Schistosoma mansoni
        "7209",  # Loa loa
        "121225",  # Pediculus humanus (body louse; typhus/trench fever vector)
        # --- Fungi ---
        "4932",  # Saccharomyces cerevisiae
        "4896",  # Schizosaccharomyces pombe
        "5476",  # Candida albicans
        "5141",  # Neurospora crassa
        "5207",  # Cryptococcus neoformans
        "746128",  # Aspergillus fumigatus
        "33169",  # Eremothecium gossypii
        "209285",  # Thermochaetoides thermophila
        # --- Plants ---
        "3702",  # Arabidopsis thaliana
        "4530",  # Oryza sativa
        "4577",  # Zea mays
        "3847",  # Glycine max
        "4081",  # Solanum lycopersicum
        "4565",  # Triticum aestivum
        "3055",  # Chlamydomonas reinhardtii
        # --- Protists ---
        "44689",  # Dictyostelium discoideum
        "5833",  # Plasmodium falciparum
        "5691",  # Trypanosoma brucei
        "5664",  # Leishmania major
        "5811",  # Toxoplasma gondii
        "5741",  # Giardia duodenalis
        "5722",  # Trichomonas vaginalis
        "5759",  # Entamoeba histolytica
        "5911",  # Tetrahymena thermophila (model ciliate; telomere/telomerase biology)
        # --- Bacteria: model organisms, pathogens, and gut microbiome ---
        "562",  # Escherichia coli
        "1423",  # Bacillus subtilis
        "1280",  # Staphylococcus aureus
        "1773",  # Mycobacterium tuberculosis
        "287",  # Pseudomonas aeruginosa
        "28901",  # Salmonella enterica
        "210",  # Helicobacter pylori
        "1313",  # Streptococcus pneumoniae
        "1639",  # Listeria monocytogenes
        "1496",  # Clostridioides difficile
        "817",  # Bacteroides fragilis
        "818",  # Bacteroides thetaiotaomicron
        "216816",  # Bifidobacterium longum
        "1590",  # Lactiplantibacillus plantarum
        "239935",  # Akkermansia muciniphila
        "853",  # Faecalibacterium prausnitzii
        "573",  # Klebsiella pneumoniae
        "1351",  # Enterococcus faecalis
        "1352",  # Enterococcus faecium
        "1314",  # Streptococcus pyogenes
        "666",  # Vibrio cholerae
        "485",  # Neisseria gonorrhoeae
        "139",  # Borreliella burgdorferi
        "813",  # Chlamydia trachomatis
        "197",  # Campylobacter jejuni
        "1502",  # Clostridium perfringens
        "1491",  # Clostridium botulinum
        "1392",  # Bacillus anthracis
        "160",  # Treponema pallidum
        "2104",  # Mycoplasmoides pneumoniae
        "235",  # Brucella abortus
        "470",  # Acinetobacter baumannii
        "727",  # Haemophilus influenzae
        "1747",  # Cutibacterium acnes
        "1360",  # Lactococcus lactis
        "274",  # Thermus thermophilus
        "624",  # Shigella sonnei
        "623",  # Shigella flexneri
        "263",  # Francisella tularensis
        "1772",  # Mycolicibacterium smegmatis (standard non-pathogenic M. tuberculosis surrogate)
        # --- Archaea ---
        "2190",  # Methanocaldococcus jannaschii
        "2242",  # Halobacterium salinarum
        "2246",  # Haloferax volcanii
        "311400",  # Thermococcus kodakarensis
        "2261",  # Pyrococcus furiosus
        # --- Viruses of human health interest ---
        "11676",  # Human immunodeficiency virus 1
        "2697049",  # SARS-CoV-2
        "11320",  # Influenza A virus
        "10407",  # Hepatitis B virus
        "10376",  # Human gammaherpesvirus 4 (Epstein-Barr)
        "10359",  # Human betaherpesvirus 5 (cytomegalovirus)
        "10298",  # Human alphaherpesvirus 1 (HSV-1)
        "333760",  # Human papillomavirus 16
        "11723",  # Simian immunodeficiency virus
        "3431487",  # Orthopoxvirus variola
        "10245",  # Orthopoxvirus vaccinia
        "3050293",  # Simplexvirus humanalpha2 (HSV-2)
        "3431483",  # Orthopoxvirus monkeypox
        "10710",  # Lambdavirus lambda (phage lambda; foundational molecular-biology model)
    }
)
