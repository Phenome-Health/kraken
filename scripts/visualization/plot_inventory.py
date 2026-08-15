#!/usr/bin/env python3
"""
Figure 2 for the KRAKEN paper: a two-panel inventory of the knowledge graph's
content, derived from a metagraph JSON file.

  Left panel  : primary knowledge sources, by edge count.
  Right panel : node identifier prefixes, by node count.

Each panel is sorted by count and shows the top N entries (default 100); the
panel header reports the true total (e.g. "showing top 100 of 106").

Three versions are written each run, differing only in the x-axis scale:
  *_linear : true proportions (small bars nearly vanish given the range).
  *_sqrt   : square-root axis -- a soft compression.
  *_log    : log axis -- strong compression so small bars stay visible.

Usage:
    python plot_inventory.py METAGRAPH.json [-o OUT.pdf] [--top N]
"""
import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import seaborn as sns

# --- single bar color (colorblind-safe, prints well) ----------------------
BAR_COLOR = "#3E7CB1"  # calm steel blue

# above this many bars, per-bar value labels become illegible at print scale,
# so they are auto-suppressed (the axis still conveys magnitude)
VALUE_LABEL_MAX = 45

# x-axis scales offered (log is the classic view; linear/sqrt for comparison)
SCALE_SUFFIX = {"linear": "", "sqrt": " (sqrt scale)", "log": " (log scale)"}


def sqrt_ticks(vmax: float, n: int = 5) -> list[float]:
    """Ticks evenly spaced in sqrt space (so they don't bunch up), each rounded
    to one significant figure for clean labels."""
    if vmax <= 0:
        return [0]
    smax = math.sqrt(vmax)
    ticks = [0]
    for i in range(1, n + 1):
        raw = (smax * i / n) ** 2
        mag = 10 ** math.floor(math.log10(raw))
        ticks.append(round(raw / mag) * mag)
    return sorted(set(ticks))


def apply_xscale(ax, scale: str, vmax: float):
    if scale == "log":
        ax.set_xscale("log")
        ax.set_xlim(left=0.8)  # log can't start at 0
    elif scale == "sqrt":
        ax.set_xscale("function", functions=(np.sqrt, np.square))
        ax.set_xlim(left=0)
        ax.set_xticks(sqrt_ticks(vmax))  # avoid the labels bunching at the high end
    else:  # linear
        ax.set_xscale("linear")
        ax.set_xlim(left=0)


def clean_label(key: str) -> str:
    """Strip the 'infores:' scheme for readability; leave prefixes as-is."""
    return key[len("infores:") :] if key.startswith("infores:") else key


def draw_panel(
    ax, data: dict, name: str, base_xlabel: str, top: int | None, value_labels: bool, panel_label: str, scale: str
):
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
    bars = ax.barh(list(y), values, color=BAR_COLOR, height=0.78, edgecolor="white", linewidth=0.3)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()  # largest at top; each panel fills its own height

    apply_xscale(ax, scale, max(values) if values else 1)
    ax.set_xlabel(base_xlabel + SCALE_SUFFIX[scale], fontsize=9)
    ax.set_title(title, fontsize=10, loc="left", pad=6)
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: "0" if v == 0 else (f"{int(v):,}" if v >= 1 else ""))
    )
    ax.tick_params(axis="x", labelsize=7.5)
    ax.grid(axis="x", which="major", linewidth=0.5, alpha=0.5)
    ax.grid(axis="x", which="minor", linewidth=0.3, alpha=0.25)
    ax.set_axisbelow(True)

    # bold panel label (A/B) in the top-left margin, for caption cross-reference
    ax.annotate(
        panel_label,
        xy=(0, 1),
        xycoords="axes fraction",
        xytext=(-38, 16),
        textcoords="offset points",
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    # value labels auto-suppress when too dense to stay legible at print scale.
    # bar_label pads in points, so it sits just past the bar end on any scale.
    if value_labels and shown <= VALUE_LABEL_MAX:
        ax.bar_label(bars, labels=[f"{v:,}" for v in values], padding=3, fontsize=6, color="#334155")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metagraph", type=Path, help="Path to metagraph JSON")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("sources.pdf"),
        help="Output base path; _linear/_sqrt/_log suffixes are added",
    )
    ap.add_argument("--top", type=int, default=100, help="Show only the top N per panel (default: 100; use 0 for all)")
    ap.add_argument("--no-value-labels", action="store_true", help="Hide the numeric labels at bar ends")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()
    # write outputs next to this script (not the cwd) unless an absolute path is given
    if not args.output.is_absolute():
        args.output = Path(__file__).resolve().parent / args.output

    meta = json.loads(args.metagraph.read_text())
    pks = meta["primary_knowledge_sources"]
    prefixes = meta["node_prefixes"]

    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
    # embed TrueType (Type 42) rather than Type 3 fonts, which journals reject
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    n_rows = max(len(pks), len(prefixes))
    if args.top and args.top > 0:
        n_rows = min(args.top, n_rows)
    height = max(4.0, 0.20 * n_rows + 1.0)

    stem, suf = args.output.with_suffix(""), args.output.suffix or ".pdf"
    for scale in ("linear", "sqrt", "log"):
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, height))
        draw_panel(axL, pks, "Primary knowledge sources", "Edges", args.top, not args.no_value_labels, "A", scale)
        draw_panel(axR, prefixes, "Node identifier prefixes", "Nodes", args.top, not args.no_value_labels, "B", scale)
        fig.tight_layout()
        out = Path(f"{stem}_{scale}{suf}")
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
        png = out.with_suffix(".png")
        if png != out:
            fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out}" + (f" and {png}" if png != out else ""))


if __name__ == "__main__":
    main()
