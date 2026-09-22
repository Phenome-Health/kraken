"""Babel clique outliers: an id whose NAME says Babel put it in the wrong clique.

Babel builds some cliques out of a shared structure, and a class-level entry that carries an example structure can
land in the wrong one: HMDB's "Glycolipids" (HMDB0302365) sits in Babel's sphingomyelin clique because the two
share an InChIKey. Nothing else can pull it out -- a clique links it to every member, and the aggregators that put
it with glycolipid (kg2, ROBOKOP) are ignored where Babel knows both ids.

The name is the tell. An id is an OUTLIER when:

* its Babel clique is name-coherent: at least two of its other members share a name (sphingomyelin's are
  "Sphingomyelins" / "sphingomyelin"),
* its own name matches none of them, and
* its name matches at least two ids of ONE other Babel clique ("Glycolipids" = MeSH's and UMLS's, both in
  glycolipid's clique).

Only name matches the pairwise guardrails allow count, and only chemicals, diseases, anatomy and the like are judged
-- not organisms, genes, proteins or variants, whose names are symbols. For an outlier, Babel's evidence about it is
dropped and the aggregator claims Babel's say had overruled come back, so its name and the other sources decide.
Measured on the 2.1.1 inputs: 2,258 outliers. In a sample of 45, about 60% are fixes (ursodeoxycholic acid in
tauroursodiol's clique, dexamfetamine in racemic amphetamine's, DOID "generalized epilepsy" in MONDO "idiopathic
generalized epilepsy"), about 30% move an id into another Babel clique of the same substance, and about 10% move it
into a salt, charge form or class whose clique holds ids of the same name (NAN 190 into NAN 190 hydrobromide).

Everything here is external sorts and merge joins over temp files, one id at a time, so memory stays flat.
"""

from __future__ import annotations

import logging
import os
import subprocess
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from kraken.entity_resolution.families import ALL_FAMILIES
from kraken.entity_resolution.guardrails import GuardrailConfig, NodeInfo, cluster_violations
from kraken.entity_resolution.name_sim import SPACED_NAME_FAMILIES
from kraken.utils.kg_io import remove_file

SEP = "\t"
WILDCARD_SIGNATURE = "*"
# An outlier needs this many name matches in one other clique: one could be a coincidence.
MIN_MATCHES_IN_OTHER_CLIQUE = 2


def families_signature(branches: frozenset[str]) -> str:
    """A node's branch families as one TSV-safe field (``*`` for a wildcard)."""
    return WILDCARD_SIGNATURE if branches is ALL_FAMILIES else ",".join(sorted(branches))


def _families(signature: str) -> frozenset[str]:
    return ALL_FAMILIES if signature == WILDCARD_SIGNATURE else frozenset(signature.split(","))


def compatible_name_pairs(members: dict[str, str], config: GuardrailConfig) -> Iterator[tuple[str, str]]:
    """The pairs of one name group (``curie -> families signature``) that the pairwise guardrails allow."""
    ids = sorted(members)
    info = {c: NodeInfo(curie=c, branches=_families(members[c]), taxon=None) for c in ids}
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if not cluster_violations([a, b], {a: info[a], b: info[b]}, config):
                yield a, b


def _sort(path: Path, temp_dir: Path, *, unique_pairs: bool = False) -> Path:
    """Byte-order sort on the first field; with ``unique_pairs``, on the first two, dropping repeats of them (``-u``
    compares only the sort keys)."""
    out = path.with_suffix(".sorted")
    keys = ["-k1,1", "-k2,2", "-u"] if unique_pairs else ["-k1,1"]
    args = ["sort", "-t", SEP, *keys, "-T", str(temp_dir), "-o", str(out), str(path)]
    subprocess.run(args, check=True, env={**os.environ, "LC_ALL": "C"})
    return out


def _join_clique(rows_path: Path, cliques_path: Path) -> Iterator[tuple[list[str], str, str]]:
    """Merge-join rows (sorted on their first field) with ``id, hub, size`` (sorted on id): each row with its first
    field's clique hub and size ("" and "0" when Babel has it in no clique)."""
    with open(rows_path) as rows, open(cliques_path) as cliques:
        clique = next(cliques, None)
        for line in rows:
            fields = line.rstrip("\n").split(SEP)
            while clique is not None and clique.split(SEP, 1)[0] < fields[0]:
                clique = next(cliques, None)
            if clique is not None and clique.split(SEP, 1)[0] == fields[0]:
                _id, hub, size = clique.rstrip("\n").split(SEP)
                yield fields, hub, size
            else:
                yield fields, "", "0"


def find_babel_outliers(
    name_pairs_path: Path, cliques_path: Path, clique_cap: int, temp_dir: Path
) -> dict[str, tuple[str, str, int]]:
    """``outlier id -> (its clique hub, the other clique's hub, how many of its names match there)``.

    ``name_pairs_path`` holds guardrail-compatible name matches as ``a, b, families(a), families(b)``;
    ``cliques_path`` holds every Babel clique member as ``id, hub, clique size``. An outlier in a clique larger than
    ``clique_cap`` is skipped: such a clique is a star from its hub, and dropping an id could disconnect it.
    """
    temps: list[Path] = []

    def temp(name: str) -> Path:
        path = temp_dir / name
        temps.extend([path, path.with_suffix(".sorted")])
        return path

    try:
        directed = temp("er_s1b_outlier_directed.tmp")
        with open(name_pairs_path) as fin, open(directed, "w") as out:
            for line in fin:
                a, b, fam_a, fam_b = line.rstrip("\n").split(SEP)
                out.write(f"{a}{SEP}{b}{SEP}{fam_a}\n{b}{SEP}{a}{SEP}{fam_b}\n")
        sorted_directed = _sort(directed, temp_dir, unique_pairs=True)
        temps.append(cliques_path.with_suffix(".sorted"))
        sorted_cliques = _sort(cliques_path, temp_dir)

        # x, partner, families(x) -> partner, x, hub(x), size(x), families(x)
        step1 = temp("er_s1b_outlier_step1.tmp")
        with open(step1, "w") as out:
            for (x, partner, fam_x), hub, size in _join_clique(sorted_directed, sorted_cliques):
                out.write(f"{partner}{SEP}{x}{SEP}{hub}{SEP}{size}{SEP}{fam_x}\n")
        sorted_step1 = _sort(step1, temp_dir)

        # -> x, hub(x), size(x), families(x), hub(partner)
        step2 = temp("er_s1b_outlier_step2.tmp")
        intra: set[str] = set()  # cliques where two members' names match
        with open(step2, "w") as out:
            for (_partner, x, hub_x, size_x, fam_x), hub_p, _size_p in _join_clique(sorted_step1, sorted_cliques):
                if hub_x and hub_x == hub_p:
                    intra.add(hub_x)
                out.write(f"{x}{SEP}{hub_x}{SEP}{size_x}{SEP}{fam_x}{SEP}{hub_p}\n")
        sorted_step2 = _sort(step2, temp_dir)

        outliers: dict[str, tuple[str, str, int]] = {}
        current = None
        partner_hubs: Counter = Counter()
        hub_x = fam_x = ""
        size_x = 0

        def decide() -> None:
            if current is None or not hub_x or hub_x not in intra or size_x > clique_cap:
                return
            if _families(fam_x) is not ALL_FAMILIES and _families(fam_x) & SPACED_NAME_FAMILIES:
                return
            if partner_hubs.get(hub_x):
                return  # its name matches its own clique
            other = [(n, h) for h, n in partner_hubs.items() if h]
            if other:
                count, hub = max(other)
                if count >= MIN_MATCHES_IN_OTHER_CLIQUE:
                    outliers[current] = (hub_x, hub, count)

        with open(sorted_step2) as fin:
            for line in fin:
                x, hub, size, fam, hub_p = line.rstrip("\n").split(SEP)
                if x != current:
                    decide()
                    current, hub_x, size_x, fam_x, partner_hubs = x, hub, int(size), fam, Counter()
                partner_hubs[hub_p] += 1
            decide()
        return outliers
    finally:
        for path in temps:
            remove_file(path)


def log_babel_outliers(outliers: dict[str, tuple[str, str, int]], report_path: Path | None = None) -> None:
    logging.info(
        "entity_resolution: %d Babel clique outliers (name matches nothing in its own clique, and %d+ ids of "
        "another); their Babel evidence is dropped and the aggregator claims about them restored",
        len(outliers),
        MIN_MATCHES_IN_OTHER_CLIQUE,
    )
    if report_path is not None and outliers:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w") as out:
            out.write(f"id{SEP}babel_clique_hub{SEP}name_matches_clique_hub{SEP}name_matches\n")
            for curie, (hub, other, count) in sorted(outliers.items()):
                out.write(f"{curie}{SEP}{hub}{SEP}{other}{SEP}{count}\n")
        logging.info("entity_resolution: Babel clique outliers listed in %s", report_path)
