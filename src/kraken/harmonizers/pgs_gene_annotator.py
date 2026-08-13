# pgs_gene_annotator.py
"""Maps GRCh38 genomic positions to overlapping protein-coding genes, for PGS -> gene linking.

'Overlap' means the variant falls within the gene body -- a high-confidence "located in" relationship, NOT a
causal one (regulatory/intergenic variants that act on distal genes are not captured). This is the pragmatic
first-pass variant-to-gene mapping; causal V2G (eQTL/colocalization, Open Targets L2G) is out of scope.
"""
import bisect
import gzip
import logging
import re
from collections import defaultdict
from pathlib import Path

# Ensembl gene model (GRCh38), pinned for reproducibility. Bump deliberately (it changes gene coordinates).
ENSEMBL_GTF_URL = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"

_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
_GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
_MAX_GENE_SPAN = 3_000_000  # a gene overlapping pos starts within this many bp before it (longest gene ~2.4 Mb)


class GeneAnnotator:
    """Overlap lookup from a GRCh38 position to protein-coding gene(s), built from an Ensembl GTF."""

    def __init__(self, gtf_path: Path):
        self._records: dict[str, list[tuple[int, int, str, str | None]]] = {}  # contig -> sorted [(start,end,ensg,sym)]
        self._starts: dict[str, list[int]] = {}  # contig -> parallel list of starts (for bisect)
        self._load(gtf_path)

    def _load(self, gtf_path: Path):
        by_contig: dict[str, list[tuple[int, int, str, str | None]]] = defaultdict(list)
        n = 0
        opener = gzip.open if str(gtf_path).endswith(".gz") else open
        with opener(gtf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 9 or cols[2] != "gene" or 'gene_biotype "protein_coding"' not in cols[8]:
                    continue
                ensg = _GENE_ID_RE.search(cols[8])
                if not ensg:
                    continue
                sym = _GENE_NAME_RE.search(cols[8])
                by_contig[cols[0]].append((int(cols[3]), int(cols[4]), ensg.group(1), sym.group(1) if sym else None))
                n += 1
        for contig, recs in by_contig.items():
            recs.sort()
            self._records[contig] = recs
            self._starts[contig] = [r[0] for r in recs]
        logging.info(f"GeneAnnotator loaded {n:,} protein-coding genes across {len(by_contig)} contigs")

    def genes_at(self, chrom: str, pos: int) -> list[tuple[str, str | None]]:
        """Return [(ensembl_gene_id, symbol)] for protein-coding genes whose body contains `pos`."""
        contig = str(chrom)
        recs = self._records.get(contig)
        if not recs:
            return []
        starts = self._starts[contig]
        # genes overlapping pos have start in [pos - _MAX_GENE_SPAN, pos]; among that window keep those with end >= pos
        lo = bisect.bisect_left(starts, pos - _MAX_GENE_SPAN)
        hi = bisect.bisect_right(starts, pos)
        return [(ensg, sym) for (start, end, ensg, sym) in recs[lo:hi] if end >= pos]
