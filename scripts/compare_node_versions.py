"""Compare how the SAME entities resolved across two KRAKEN builds.

Point it at two nodes files (e.g. 2.1.0 vs a fresh build). For each watchlist id
(or --ids) it finds the node in BOTH files -- matching on representative id OR
equivalent_ids membership, so it's robust to the canonical id changing between
builds -- and shows the membership delta: how many ids merged in, which appeared
or disappeared, and whether the category set changed. This is the quick "did the
disease cliques re-form / did anything over-merge" check.

    uv run python scripts/compare_node_versions.py OLD_NODES NEW_NODES
    uv run python scripts/compare_node_versions.py \
        /Volumes/AmySSD/kraken-data/xx_artifacts_v2.1.0/integrated/kraken_nodes_2.1.0.jsonl \
        /Volumes/AmySSD/kraken-data/artifacts/integrated/kraken_nodes_2.1.1.jsonl
    uv run python scripts/compare_node_versions.py OLD NEW --ids MONDO:0005180,NCBIGene:7157
    uv run python scripts/compare_node_versions.py OLD NEW --counts        # also total node counts (slow)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from peek_nodes import DEFAULT_WATCHLIST, Style, find_nodes  # noqa: E402


def _members(node: dict | None) -> set[str]:
    if not node:
        return set()
    return {node.get("id")} | set(node.get("equivalent_ids") or [])


def _fmt_ids(ids: set[str], limit: int = 20) -> str:
    shown = sorted(ids)
    if len(shown) > limit:
        return ", ".join(shown[:limit]) + f"  …(+{len(shown) - limit} more)"
    return ", ".join(shown)


def _count_lines(path: Path) -> int:
    n = 0
    with open(path, "rb") as fh:
        for _ in fh:
            n += 1
    return n


def compare_one(label: str, wid: str, old: dict | None, new: dict | None, s: Style, width: int) -> None:
    print(s.dim("─" * width))
    print(f"{s.green('●')} {s.bold(label)}   {s.dim(f'[watch: {wid}]')}")

    if old is None and new is None:
        print(f"  {s.red('absent in BOTH builds')}")
        return
    if old is None:
        print(f"  {s.yellow('NEW ONLY')} — absent in old build")
    if new is None:
        print(f"  {s.red('GONE')} — present in old build, absent in new")

    old_m, new_m = _members(old), _members(new)
    old_cats = sorted(old.get("categories") or []) if old else []
    new_cats = sorted(new.get("categories") or []) if new else []

    def line(tag: str, node: dict | None, members: set[str], cats: list[str]) -> None:
        if node is None:
            print(f"  {s.dim(tag)}  {s.dim('—')}")
            return
        print(
            f"  {s.dim(tag)}  rep={s.cyan(node.get('id', '?'))}  "
            f"members={s.bold(str(len(members)))}  cats=[{s.green(', '.join(cats))}]"
        )

    line("old", old, old_m, old_cats)
    line("new", new, new_m, new_cats)

    if old and new:
        delta = len(new_m) - len(old_m)
        dcolor = s.green if delta == 0 else (s.yellow if delta > 0 else s.red)
        print(f"  {s.dim('Δ')}    members {dcolor(f'{delta:+d}')}"
              + (f"   {s.dim('(rep changed)')}" if old.get("id") != new.get("id") else ""))

        added = new_m - old_m
        removed = old_m - new_m
        if added:
            print(f"  {s.yellow('+')}    gained in new ({len(added)}): {s.grey(_fmt_ids(added))}")
        if removed:
            print(f"  {s.red('-')}    lost in new ({len(removed)}): {s.grey(_fmt_ids(removed))}")
        if not added and not removed:
            print(f"  {s.green('=')}    identical membership")
        if old_cats != new_cats:
            print(f"  {s.magenta('!')}    categories changed: {s.red(', '.join(old_cats) or '—')} "
                  f"{s.dim('→')} {s.green(', '.join(new_cats) or '—')}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("old", type=Path, help="OLD nodes .jsonl")
    ap.add_argument("new", type=Path, help="NEW nodes .jsonl")
    ap.add_argument("--ids", default=None, help="comma-separated curies (default: the watchlist)")
    ap.add_argument("--counts", action="store_true", help="also print total node counts (streams both files fully)")
    ap.add_argument("--no-color", action="store_true")
    args = ap.parse_args()

    s = Style(enabled=not args.no_color)
    width = 110
    for p in (args.old, args.new):
        if not p.exists():
            raise SystemExit(f"nodes file not found: {p}")

    watchlist = (
        [(cid, cid) for cid in (x.strip() for x in args.ids.split(",")) if cid] if args.ids else DEFAULT_WATCHLIST
    )
    watch_ids = {cid for _l, cid in watchlist}

    print(s.dim(f"OLD {args.old}"))
    print(s.dim(f"NEW {args.new}"))
    old_found = find_nodes(args.old, watch_ids)
    new_found = find_nodes(args.new, watch_ids)

    for label, cid in watchlist:
        compare_one(label, cid, old_found.get(cid), new_found.get(cid), s, width)
    print(s.dim("─" * width))

    if args.counts:
        print(s.dim("counting total nodes (streaming both files)…"))
        old_n, new_n = _count_lines(args.old), _count_lines(args.new)
        delta = new_n - old_n
        dcolor = s.yellow if delta > 0 else (s.red if delta < 0 else s.green)
        print(f"{s.bold('Total nodes')}:  old={old_n:,}   new={new_n:,}   Δ={dcolor(f'{delta:+,}')}")


if __name__ == "__main__":
    main()
