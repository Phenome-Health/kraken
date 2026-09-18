"""Audit src/kraken/harmonizers/ncbigene_taxon_allowlist.py against a built KRAKEN node set.

The allowlist scopes the NCBI Gene ingest (NCBI's all-species `gene_info` is ~72M genes across ~54k taxa).
It is deliberately CURATED, not derived -- see that module's docstring for the measurements behind that
choice. This script does not rewrite it. It answers two questions:

  1. What fraction of the NCBIGene ids already in a build does the curated list cover?
  2. Which uncovered species account for the most of them -- i.e. what is worth ADDING?

Deriving the list from those answers would be circular (it could never introduce an organism the graph
doesn't already have). Use them as evidence, then decide which organisms genuinely belong.

    uv run python scripts/audit_ncbigene_taxon_allowlist.py \
        --kraken-nodes /path/to/kraken_nodes_2.1.0.jsonl \
        --gene-info /path/to/input_data/ncbigene/gene_info.gz \
        --taxdump /path/to/input_data/ncbigene/taxdump.tar.gz
"""

import argparse
import gzip
import logging
import re
from collections import Counter
from pathlib import Path

from kraken.harmonizers.helpers.ncbigene_taxon_allowlist import TAXON_ALLOWLIST
from kraken.utils.taxonomy import TaxonNormalizer

# Matched against raw JSONL text rather than parsed JSON: an NCBIGene CURIE is equally interesting whether it
# appears as a node's `id` or inside its `equivalent_ids`, and scanning text avoids parsing ~8GB of JSON.
NCBIGENE_CURIE_PATTERN = re.compile(r'"NCBIGene:(\d+)"')

GENE_INFO_TAX_ID_COLUMN = 0
GENE_INFO_GENE_ID_COLUMN = 1

UNCOVERED_TO_REPORT = 25


def _open_maybe_gzipped(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def collect_gene_ids(kraken_nodes_path: Path) -> set[str]:
    """Every distinct NCBI GeneID referenced by a KRAKEN nodes file, from `id` or `equivalent_ids` alike."""
    logging.info(f"Scanning {kraken_nodes_path} for NCBIGene CURIEs..")
    gene_ids: set[str] = set()
    with _open_maybe_gzipped(kraken_nodes_path) as nodes_file:
        for line_num, line in enumerate(nodes_file, 1):
            if line_num % 1_000_000 == 0:
                logging.info(f"    at {line_num} nodes ({len(gene_ids)} gene ids so far)")
            gene_ids.update(NCBIGENE_CURIE_PATTERN.findall(line))
    logging.info(f"Found {len(gene_ids)} distinct NCBIGene ids")
    return gene_ids


def collect_species(gene_info_path: Path, gene_ids: set[str], taxonomy: TaxonNormalizer) -> tuple[Counter, Counter]:
    """Return (genes-per-species among the build's gene ids, total genes per species in gene_info). The
    second is the ingest cost of allowlisting a species."""
    logging.info(f"Looking up taxa in {gene_info_path}..")
    build_species: Counter = Counter()
    all_species: Counter = Counter()
    with _open_maybe_gzipped(gene_info_path) as gene_info_file:
        gene_info_file.readline()  # header
        for line_num, line in enumerate(gene_info_file, 1):
            if line_num % 10_000_000 == 0:
                logging.info(f"    at {line_num} gene_info rows")
            columns = line.split("\t", GENE_INFO_GENE_ID_COLUMN + 1)
            if len(columns) <= GENE_INFO_GENE_ID_COLUMN:
                continue
            species = taxonomy.to_species(columns[GENE_INFO_TAX_ID_COLUMN])
            all_species[species] += 1
            if columns[GENE_INFO_GENE_ID_COLUMN] in gene_ids:
                build_species[species] += 1
    return build_species, all_species


def report(build_species: Counter, all_species: Counter, taxonomy: TaxonNormalizer) -> None:
    allowed = {taxonomy.to_species(tax_id) for tax_id in TAXON_ALLOWLIST}
    total = sum(build_species.values())
    covered = sum(count for species, count in build_species.items() if species in allowed)
    ingest = sum(all_species[species] for species in allowed)

    print()
    print(f"Allowlist: {len(TAXON_ALLOWLIST)} curated organisms -> {len(allowed)} species after rollup")
    print(f"Coverage:  {covered:,} of {total:,} of the build's NCBIGene ids ({100 * covered / total:.1f}%)")
    print(f"Ingest:    {ingest:,} genes from gene_info")

    uncovered = [(s, c) for s, c in build_species.most_common() if s not in allowed]
    if not uncovered:
        print("\nEvery species in the build is covered.")
        return
    print(f"\nBiggest uncovered species ({len(uncovered)} total) -- candidates to ADD if they belong:")
    for species, count in uncovered[:UNCOVERED_TO_REPORT]:
        name = taxonomy.scientific_name(species) or "?"
        print(f"  {species:>9}  {name:<45} {count:>6} genes in build, {all_species[species]:>7} to ingest")


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kraken-nodes", type=Path, required=True, help="A built KRAKEN nodes .jsonl (or .gz)")
    parser.add_argument("--gene-info", type=Path, required=True, help="NCBI's gene_info file (or .gz)")
    parser.add_argument("--taxdump", type=Path, required=True, help="NCBI's taxdump.tar.gz")
    args = parser.parse_args()

    gene_ids = collect_gene_ids(args.kraken_nodes)
    if not gene_ids:
        raise SystemExit(f"No NCBIGene CURIEs found in {args.kraken_nodes} - is that a KRAKEN nodes file?")

    taxonomy = TaxonNormalizer(args.taxdump)
    build_species, all_species = collect_species(args.gene_info, gene_ids, taxonomy)
    if not build_species:
        raise SystemExit(f"None of the {len(gene_ids)} gene ids appear in {args.gene_info}")

    report(build_species, all_species, taxonomy)


if __name__ == "__main__":
    main()
