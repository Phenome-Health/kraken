#!/usr/bin/env python3
"""
Figure 1 for the KRAKEN paper: a two-panel inventory of the knowledge graph's
content, derived from a metagraph JSON file.

  Left panel  : primary knowledge sources, by edge count (log axis).
  Right panel : node identifier prefixes, by node count (log axis).

Each panel is sorted by count and shows the top N entries (default 100); the
panel header reports the true total (e.g. "showing top 100 of 106").

Usage:
    python plot_inventory.py METAGRAPH.json [-o OUT.pdf] [--top N]
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

# --- single bar color (colorblind-safe, prints well) ----------------------
BAR_COLOR = "#3E7CB1"   # calm steel blue

# above this many bars, per-bar value labels become illegible at print scale,
# so they are auto-suppressed (the log axis still conveys magnitude)
VALUE_LABEL_MAX = 45


def clean_label(key: str) -> str:
    """Strip the 'infores:' scheme for readability; leave prefixes as-is."""
    return key[len("infores:"):] if key.startswith("infores:") else key


def draw_panel(ax, data: dict, name: str, xlabel: str,
               top: int | None, value_labels: bool, panel_label: str):
    total = len(data)
    items = sorted(data.items(), key=lambda kv: kv[1], reverse=True)
    if top and top > 0:
        items = items[:top]
    shown = len(items)
    if shown < total:
        title = f"{name} (showing top {shown} of {total:,})"
    else:
        title = f"{name} (n={total:,})"
    labels = [clean_label(k) for k, _ in items]
    values = [v for _, v in items]

    y = range(len(items))
    ax.barh(list(y), values, color=BAR_COLOR, height=0.78,
            edgecolor="white", linewidth=0.3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()  # largest at top; each panel fills its own height

    ax.set_xscale("log")
    ax.set_xlim(left=0.8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=10, loc="left", pad=6)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{int(v):,}" if v >= 1 else ""))
    ax.tick_params(axis="x", labelsize=7.5)
    ax.grid(axis="x", which="major", linewidth=0.5, alpha=0.5)
    ax.grid(axis="x", which="minor", linewidth=0.3, alpha=0.25)
    ax.set_axisbelow(True)

    # bold panel label (A/B) in the top-left margin, for caption cross-reference
    ax.annotate(panel_label, xy=(0, 1), xycoords="axes fraction",
                xytext=(-38, 16), textcoords="offset points",
                fontsize=14, fontweight="bold", ha="left", va="bottom")

    # value labels auto-suppress when too dense to stay legible at print scale
    if value_labels and shown <= VALUE_LABEL_MAX:
        for yi, v in zip(y, values):
            ax.text(v * 1.08, yi, f"{v:,}", va="center", ha="left",
                    fontsize=6, color="#334155")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metagraph", type=Path, help="Path to metagraph JSON")
    ap.add_argument("-o", "--output", type=Path, default=Path("figure1_inventory.pdf"),
                    help="Output figure path (extension sets format; default PDF)")
    ap.add_argument("--top", type=int, default=100,
                    help="Show only the top N per panel (default: 100; use 0 for all)")
    ap.add_argument("--no-value-labels", action="store_true",
                    help="Hide the numeric labels at bar ends")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    meta = json.loads(args.metagraph.read_text())
    pks = meta["primary_knowledge_sources"]
    prefixes = meta["node_prefixes"]

    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams["font.family"] = "sans-serif"
    # embed TrueType (Type 42) rather than Type 3 fonts, which journals reject
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    n_rows = max(len(pks), len(prefixes))
    if args.top and args.top > 0:
        n_rows = min(args.top, n_rows)
    height = max(4.0, 0.20 * n_rows + 1.0)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, height))

    draw_panel(axL, pks, "Primary knowledge sources",
               "Edges (log scale)", args.top, not args.no_value_labels, "A")
    draw_panel(axR, prefixes, "Node identifier prefixes",
               "Nodes (log scale)", args.top, not args.no_value_labels, "B")

    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight")
    png = args.output.with_suffix(".png")
    if png != args.output:
        fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
    print(f"Wrote {args.output}" + (f" and {png}" if png != args.output else ""))


if __name__ == "__main__":
    main()
