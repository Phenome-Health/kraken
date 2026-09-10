"""Pretty-print a handful of well-known "watchlist" nodes from a KRAKEN nodes file.

A quick eyeball check after a build: pull a few nodes you know well (Parkinson's,
type 2 diabetes, TP53, ...) by a preferred-vocabulary id and print every property
in a scannable terminal layout, so you can glance and confirm nothing looks majorly
wrong (did the disease clique merge? are the categories sane? did the gene/protein
conflate?).

A node matches a watch id if that id is its representative ``id`` OR appears in its
``equivalent_ids`` -- so it's found regardless of which member became canonical.

    uv run python scripts/peek_nodes.py                     # newest build (from build_config.yaml)
    uv run python scripts/peek_nodes.py --nodes /path/to/kraken_nodes_2.1.0.jsonl
    uv run python scripts/peek_nodes.py --ids MONDO:0005180,NCBIGene:7157
    uv run python scripts/peek_nodes.py --no-color | less -R
"""

from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from pathlib import Path

# (label, preferred-vocabulary curie). Edit freely -- this is just a spot-check set.
DEFAULT_WATCHLIST: list[tuple[str, str]] = [
    ("Parkinson disease", "MONDO:0005180"),
    ("Atrial fibrillation", "MONDO:0004981"),
    ("Type 2 diabetes mellitus", "MONDO:0005148"),
    ("Diabetes mellitus", "MONDO:0005015"),
    ("Alzheimer disease", "MONDO:0004975"),
    ("TP53 (gene/protein)", "NCBIGene:7157"),
    ("BRCA1 (gene/protein)", "NCBIGene:672"),
    ("Metformin (drug)", "CHEBI:6801"),
    ("Glucose (metabolite)", "CHEBI:17234"),
    ("Heart (anatomy)", "UBERON:0000948"),
    ("Apoptosis (process)", "GO:0006915"),
    ("Seizure (phenotype)", "HP:0001250"),
]

# Fields printed first, in this order; any remaining keys are printed after.
FIELD_ORDER = [
    "id",
    "categories",
    "name",
    "taxon",
    "provided_by",
    "description",
    "synonyms",
    "equivalent_ids",
    "urls",
    "attributes",
]


class Style:
    def __init__(self, enabled: bool):
        self.on = enabled

    def __call__(self, code: str, text: str) -> str:
        return f"\033[{code}m{text}\033[0m" if self.on else text

    def bold(self, t): return self("1", t)
    def dim(self, t): return self("2", t)
    def cyan(self, t): return self("36", t)
    def green(self, t): return self("32", t)
    def yellow(self, t): return self("33", t)
    def magenta(self, t): return self("35", t)
    def red(self, t): return self("31", t)
    def grey(self, t): return self("90", t)


def default_nodes_path() -> Path | None:
    """Derive the current build's nodes file from build_config.yaml (base_path +
    integration output dir + version), so the default tracks the latest build."""
    cfg_path = Path(__file__).resolve().parents[1] / "config" / "build_config.yaml"
    if not cfg_path.exists():
        return None
    import yaml

    cfg = yaml.safe_load(cfg_path.read_text())
    base = cfg.get("base_path") or Path(__file__).resolve().parents[1]
    out_dir = (cfg.get("integration") or {}).get("output_directory", "artifacts/integrated/")
    version = cfg.get("kraken_version", "")
    return Path(base) / out_dir / f"kraken_nodes_{version}.jsonl"


def find_nodes(nodes_path: Path, watch_ids: set[str]) -> dict[str, dict]:
    """One streaming pass; return {watch_id: node} for the first node matching each
    watch id (by representative id or equivalent_ids membership). Early-exits once
    every watch id is accounted for."""
    remaining = set(watch_ids)
    found: dict[str, dict] = {}
    with open(nodes_path) as fh:
        for line in fh:
            if not remaining:
                break
            # cheap prefilter before json.loads on a multi-GB file
            if not any(wid in line for wid in remaining):
                continue
            node = json.loads(line)
            ids = {node.get("id")} | set(node.get("equivalent_ids") or [])
            for wid in list(remaining):
                if wid in ids:
                    found[wid] = node
                    remaining.discard(wid)
    return found


def _wrap(text: str, width: int) -> list[str]:
    import textwrap

    return textwrap.wrap(text, width=width) or [""]


def _prefix_summary(ids: list[str], s: Style, top: int = 12) -> str:
    counts = Counter(i.split(":", 1)[0] for i in ids)
    parts = [f"{s.yellow(p)}:{n}" for p, n in counts.most_common(top)]
    extra = len(counts) - top
    if extra > 0:
        parts.append(s.dim(f"(+{extra} prefixes)"))
    return "  ".join(parts)


def print_node(label: str, watch_id: str, node: dict, s: Style, width: int) -> None:
    keycol = 13
    valwidth = max(20, width - keycol - 2)
    rep = node.get("id", "?")
    cats = ", ".join(node.get("categories") or []) or s.dim("(none)")

    print(s.dim("─" * width))
    header = f"{s.green('●')} {s.bold(label)}"
    tag = s.dim(f"[watch: {watch_id}]")
    pad = max(1, width - len(label) - len(watch_id) - 12)
    print(f"{header}{' ' * pad}{tag}")
    print(f"  {s.cyan(s.bold(rep))}   {s.green(cats)}")

    def row(key: str, lines: list[str], color=None) -> None:
        if not lines:
            return
        first = lines[0]
        klabel = s.dim(f"{key:<{keycol}}")
        val0 = color(first) if color else first
        print(f"  {klabel}{val0}")
        for cont in lines[1:]:
            print(f"  {' ' * keycol}{color(cont) if color else cont}")

    printed = {"id", "categories"}

    def emit(key: str) -> None:
        if key in printed or key not in node:
            return
        printed.add(key)
        val = node[key]
        if key == "equivalent_ids":
            ids = list(val)
            row(f"equiv ({len(ids)})", [_prefix_summary(ids, s)])
            joined = ", ".join(ids)
            for ln in _wrap(joined, valwidth):
                print(f"  {' ' * keycol}{s.grey(ln)}")
        elif key == "synonyms":
            syn = list(val)
            shown = ", ".join(syn[:20])
            suffix = s.dim(f"  (+{len(syn) - 20} more)") if len(syn) > 20 else ""
            lines = _wrap(shown, valwidth)
            row(f"synonyms({len(syn)})", lines)
            if suffix:
                print(f"  {' ' * keycol}{suffix}")
        elif key == "attributes":
            row("attributes", [""])
            for src, attrs in val.items():
                compact = json.dumps(attrs, ensure_ascii=False)
                if len(compact) > 160:  # glance tool: keep attribute blobs short
                    compact = compact[:160].rstrip() + s.dim(f" …(+{len(compact) - 160} chars)")
                for i, ln in enumerate(_wrap(f"{src}: {compact}", valwidth)):
                    prefix = s.magenta(ln) if i == 0 else s.dim(ln)
                    print(f"  {' ' * keycol}{prefix}")
        elif isinstance(val, list):
            row(key, _wrap(", ".join(map(str, val)), valwidth))
        elif key == "description":
            text = str(val)
            if len(text) > 280:  # descriptions concatenate across sources; trim for glancing
                text = text[:280].rstrip() + f" …(+{len(text) - 280} chars)"
            row(key, _wrap(text, valwidth), color=s.dim)
        elif key == "name":
            row(key, [s.bold(str(val))])
        else:
            row(key, _wrap(str(val), valwidth))

    for key in FIELD_ORDER:
        emit(key)
    for key in node:  # anything not in FIELD_ORDER
        emit(key)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nodes", type=Path, default=None, help="nodes .jsonl (default: current build)")
    ap.add_argument("--ids", default=None, help="comma-separated curies to peek instead of the default watchlist")
    ap.add_argument("--no-color", action="store_true")
    args = ap.parse_args()

    s = Style(enabled=not args.no_color)
    width = min(shutil.get_terminal_size((110, 24)).columns, 130)

    nodes_path = args.nodes or default_nodes_path()
    if not nodes_path or not Path(nodes_path).exists():
        raise SystemExit(f"nodes file not found: {nodes_path} (pass --nodes)")

    if args.ids:
        watchlist = [(cid, cid) for cid in (x.strip() for x in args.ids.split(",")) if cid]
    else:
        watchlist = DEFAULT_WATCHLIST

    watch_ids = {cid for _label, cid in watchlist}
    print(s.dim(f"reading {nodes_path}"))
    found = find_nodes(Path(nodes_path), watch_ids)

    for label, cid in watchlist:
        node = found.get(cid)
        if node is None:
            print(s.dim("─" * width))
            tag = s.dim(f"[watch: {cid}]")
            print(f"{s.red('✗')} {s.bold(label)}   {tag}   {s.red('NOT FOUND (absent or dropped)')}")
            continue
        print_node(label, cid, node, s, width)
    print(s.dim("─" * width))

    missing = [cid for _l, cid in watchlist if cid not in found]
    tail = f" Missing: {', '.join(missing)}" if missing else ""
    print(f"\n{s.bold('Summary')}: {len(found)}/{len(watchlist)} found.{tail}")


if __name__ == "__main__":
    main()
