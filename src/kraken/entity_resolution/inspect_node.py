"""Print everything about the node an id ended up in: every member id, what it is, and why it's there.

    uv run python -m kraken.entity_resolution.inspect_node MONDO:0005148
    uv run python -m kraken.entity_resolution.inspect_node HMDB:HMDB0302365 UMLS:C0011847
    uv run python -m kraken.entity_resolution.inspect_node CHEBI:6801 --db path/to/er_debug_2.2.0.sqlite

Any member id finds the node. For each member: Babel's name and type for it, what every other source called it, its
taxon, which sources provided it, and which Babel clique it's in. Then the evidence: what ties together the Babel
cliques inside the node (a node of one clique needs no explaining), Babel clique members that ended up in OTHER nodes
(a split), and the strongest links from members to ids outside it (a near miss, or what was cut).

Reads the ER debug database a build writes to its integrated debug dir (see ``debug_db``); by default the one for
the build configured in ``config/build_config.yaml``.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import yaml

from kraken.entity_resolution.debug_db import DebugDb, IdInfo, debug_db_path
from kraken.utils.constants import PROJECT_ROOT

OUTSIDE_LINKS_SHOWN = 15
BRIDGES_SHOWN = 25


def default_db_path() -> Path:
    cfg = yaml.safe_load((PROJECT_ROOT / "config" / "build_config.yaml").read_text())
    base = Path(cfg.get("base_path") or PROJECT_ROOT)
    out_dir = (cfg.get("integration") or {}).get("output_directory", "artifacts/integrated/")
    return debug_db_path(base / out_dir / "debug", cfg.get("kraken_version"))


def _short(category: str) -> str:
    return category.removeprefix("biolink:")


def _clique_letters(members: list[IdInfo]) -> dict[str, str]:
    """A letter per Babel clique in the node, biggest share first."""
    counts = Counter(m.babel_hub for m in members if m.babel_hub)
    letters = {}
    for i, (hub, _n) in enumerate(counts.most_common()):
        letters[hub] = chr(ord("A") + i) if i < 26 else f"#{i + 1}"
    return letters


def render(db: DebugDb, curie: str, out=sys.stdout) -> None:
    def line(text: str = "") -> None:
        print(text, file=out)

    node = db.node_of(curie)
    if node is None:
        line(f"{curie}: not in this build")
        return
    members = db.members(node)
    ids = {m.id for m in members}
    letters = _clique_letters(members)
    categories = Counter(c for m in members for c in m.categories)
    line("=" * 110)
    line(f"{node}  --  {len(members)} ids  (looked up by {curie})")
    line(f"categories: {', '.join(f'{_short(c)} x{n}' for c, n in categories.most_common())}")
    line("=" * 110)

    width = max(len(m.id) for m in members) + 2
    for m in members:
        clique = f"[{letters[m.babel_hub]}]" if m.babel_hub else "[-]"
        types = ",".join(_short(c) for c in (m.babel_categories or m.categories))
        name = m.babel_name or "(no Babel name)"
        marker = " <" if m.id == curie else ""
        line(f"{clique:<5}{m.id:<{width}}{name[:60]:<62}{types[:30]}{marker}")
        details = []
        if m.taxon:
            details.append(m.taxon)
        if m.provided_by:
            details.append("from " + ", ".join(p.removeprefix("infores:") for p in m.provided_by))
        if details:
            line(f"{'':<5}{'':<{width}}{'; '.join(details)[:100]}")
        for source, other_name in m.source_names:
            if other_name != m.babel_name:
                line(f"{'':<5}{'':<{width}}{source}: {other_name[:80]}")

    if letters:
        line()
        line("Babel cliques:")
        for hub, letter in letters.items():
            clique = db.babel_clique(hub)
            here = sum(1 for c in clique if c in ids)
            elsewhere = [c for c in clique if c not in ids]
            spread = f", {len(elsewhere)} elsewhere" if elsewhere else ""
            line(f"  [{letter}] {hub}: {len(clique)} ids, {here} here{spread}")
            for other in elsewhere[:10]:
                line(f"        {other:<34} -> node {db.node_of(other)}")
    unclustered = [m.id for m in members if not m.babel_hub]
    if unclustered and len(members) > 1:
        line(f"  [-] in no Babel clique: {len(unclustered)} ids")

    evidence = list(db.evidence_touching(ids))
    hub_of = {m.id: m.babel_hub for m in members}
    inside = [e for e in evidence if e.a in ids and e.b in ids]
    bridges = [e for e in inside if not hub_of.get(e.a) or hub_of.get(e.a) != hub_of.get(e.b)]
    line()
    line(f"Evidence between ids of different Babel cliques inside the node ({len(bridges)} pairs):")
    for e in sorted(bridges, key=lambda e: -e.weight)[:BRIDGES_SHOWN]:
        line(
            f"  {e.a} [{letters.get(hub_of.get(e.a), '-')}] -- {e.b} [{letters.get(hub_of.get(e.b), '-')}]"
            f"  {e.weight:g}{'' if e.above_tau else ' (below tau)'}  {e.kinds}"
        )
    if len(bridges) > BRIDGES_SHOWN:
        line(f"  ... {len(bridges) - BRIDGES_SHOWN} more")

    outside = [e for e in evidence if (e.a in ids) != (e.b in ids)]
    line()
    line(f"Strongest links to ids outside the node ({len(outside)} in all):")
    for e in sorted(outside, key=lambda e: -e.weight)[:OUTSIDE_LINKS_SHOWN]:
        member, other = (e.a, e.b) if e.a in ids else (e.b, e.a)
        other_node = db.node_of(other)
        line(
            f"  {member} -- {other} (node {other_node})  {e.weight:g}{'' if e.above_tau else ' (below tau)'}  {e.kinds}"
        )


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Print the node an id ended up in, with every member id's details.")
    parser.add_argument("curies", nargs="+", help="any member id of the node(s) to print")
    parser.add_argument("--db", type=Path, help="the ER debug database (default: the configured build's)")
    args = parser.parse_args(argv)
    db = DebugDb(args.db or default_db_path())
    for curie in args.curies:
        render(db, curie)
        print()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
