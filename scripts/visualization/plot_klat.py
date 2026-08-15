#!/usr/bin/env python3
"""
Provenance-composition figures for the KRAKEN paper ("Content & Scale" section):
the graph's edges characterized by knowledge level and agent type (KLAT),
derived from a metagraph JSON file.

Two views (choose with --mode; default writes both):
  bars    : two horizontal-bar panels (knowledge level, agent type) sharing one
            x-axis, sorted by edge count, labeled with count + percent, and
            colored on a single-hue scale by KLAT confidence weight (darker =
            higher confidence). Uses only the marginal distributions.
  heatmap : knowledge level x agent type as a bubble grid -- disc SIZE = edge
            count (in three scales, linear/sqrt/log), disc COLOR = the cell's
            KLAT confidence (KL x AT weight). Rows/columns ordered so the
            high-confidence corner is top-left. Requires the metagraph's
            `klat_joint` field (added to the metagraph build; regenerate an
            older metagraph).

The KLAT weights mirror the values KESTREL uses for confidence-aware ranking.

Usage:
    python plot_klat.py METAGRAPH.json [-o OUT.pdf] [--mode bars|heatmap|both] [--dpi 300]
"""
import argparse
import json
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, Normalize

# KLAT confidence weights (mirror the values KESTREL uses for ranking); higher
# = more trustworthy provenance. Unlisted values fall back to 0.5.
KL_WEIGHTS = {
    "knowledge_assertion": 1.0,
    "logical_entailment": 0.9,
    "statistical_association": 0.7,
    "observation": 0.65,
    "prediction": 0.55,
    "not_provided": 0.5,
    "text_co_occurrence": 0.4,
}
AT_WEIGHTS = {
    "manual_agent": 1.0,
    "manual_validation_of_automated_agent": 0.95,
    "computational_model": 0.65,
    "data_analysis_pipeline": 0.6,
    "automated_agent": 0.55,
    "not_provided": 0.5,
    "text_mining_agent": 0.4,
    "image_processing": 0.4,
}

# single-hue confidence colormap (light -> dark = low -> high KLAT weight); used
# by the bars (per-dimension weight) and the heatmap discs (KL x AT product)
CONF_CMAP = LinearSegmentedColormap.from_list("conf", plt.cm.Blues(np.linspace(0.30, 1.0, 256)))
CONF_NORM = Normalize(0.4, 1.0)


def human(n: float) -> str:
    """Compact count, e.g. 62_700_000 -> '62.7M'."""
    if n >= 1e6:
        return f"{n / 1e6:.1f}M"
    if n >= 1e3:
        return f"{n / 1e3:.0f}K"
    return str(int(n))


def _draw_bar_panel(ax, data: dict, weights: dict, title: str, total: int, xmax: float):
    segs = sorted(data.items(), key=lambda kv: -kv[1])
    counts = [c for _, c in segs]
    colors = [CONF_CMAP(CONF_NORM(weights.get(k, 0.5))) for k, _ in segs]
    y = range(len(segs))
    ax.barh(list(y), counts, color=colors, edgecolor="white", height=0.72)
    ax.set_yticks(list(y))
    ax.set_yticklabels([k.replace("_", " ") for k, _ in segs], fontsize=8)
    ax.invert_yaxis()
    for yi, c in zip(y, counts):
        ax.text(
            c + xmax * 0.01,
            yi,
            f"{human(c)} ({100 * c / total:.1f}%)",
            va="center",
            ha="left",
            fontsize=7.5,
            color="#333333",
        )
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.set_xlim(0, xmax)
    ax.grid(axis="y", visible=False)


def make_bars(meta: dict, out: Path, dpi: int):
    total = meta["summary"]["total_edges"]
    kl, at = meta["knowledge_levels"], meta["agent_types"]
    # shared x-max: largest bar + ~18% headroom for the value labels (matplotlib
    # places nice round ticks within this; the axis need not end on one). Headroom
    # is ~scale-invariant, so this stays tight as edges grow.
    xmax = max(list(kl.values()) + list(at.values())) * 1.18

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, figsize=(10, 4.6), sharex=True, gridspec_kw={"height_ratios": [len(kl), len(at)]}
    )
    _draw_bar_panel(ax_a, kl, KL_WEIGHTS, "(A) Knowledge level", total, xmax)
    _draw_bar_panel(ax_b, at, AT_WEIGHTS, "(B) Agent type", total, xmax)
    ax_b.xaxis.set_major_formatter(lambda v, _: f"{v / 1e6:.0f}M" if v > 0 else "0")
    ax_b.set_xlabel(f"Edges (of {human(total)} total)")

    sm = cm.ScalarMappable(norm=CONF_NORM, cmap=CONF_CMAP)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=[ax_a, ax_b], fraction=0.04, pad=0.02)
    cb.set_label("KLAT confidence weight", fontsize=8)
    _save(fig, out, dpi)


def make_heatmap(meta: dict, out: Path, dpi: int, scale: str):
    joint = meta.get("klat_joint")
    if not joint:
        raise SystemExit(
            "This metagraph has no 'klat_joint'; regenerate it with the "
            "updated metagraph build to use --mode heatmap."
        )
    # rows = knowledge levels, cols = agent types, both ordered by confidence
    # weight (descending) so the high-confidence corner is top-left.
    kls = sorted(joint, key=lambda k: -KL_WEIGHTS.get(k, 0.5))
    ats = sorted({a for row in joint.values() for a in row}, key=lambda a: -AT_WEIGHTS.get(a, 0.5))
    M = np.zeros((len(kls), len(ats)))
    for i, k in enumerate(kls):
        for a, c in joint[k].items():
            M[i, ats.index(a)] = c
    nrows, ncols, vmax = len(kls), len(ats), M.max()

    # disc SIZE = edge count under the chosen scale (area ~ transform(count));
    # disc COLOR = the cell's KLAT confidence = KL_weight x AT_weight.
    tform = {"linear": lambda x: x, "sqrt": np.sqrt, "log": np.log1p}[scale]

    def area(c):
        return 24 + 1350 * tform(c) / tform(vmax)  # scatter s (points^2), with a visible floor

    cnorm = Normalize(0.16, 1.0)  # KLAT product ranges 0.4*0.4 .. 1.0*1.0
    xf, yf, sf, cf, cnt, xe, ye = [], [], [], [], [], [], []
    for i, k in enumerate(kls):
        for j, a in enumerate(ats):
            c = M[i, j]
            if c == 0:
                xe.append(j)
                ye.append(i)
                continue
            xf.append(j)
            yf.append(i)
            sf.append(area(c))
            cnt.append(c)
            cf.append(KL_WEIGHTS.get(k, 0.5) * AT_WEIGHTS.get(a, 0.5))

    fig, ax = plt.subplots(figsize=(0.6 * ncols + 3.0, 0.6 * nrows + 2.2))
    ax.scatter(xe, ye, s=16, c="#e4e4e4", edgecolors="none", zorder=1)  # empty cells
    sc = ax.scatter(xf, yf, s=sf, c=cf, cmap=CONF_CMAP, norm=cnorm, edgecolors="white", linewidths=0.6, zorder=2)
    for x, y, c, conf in zip(xf, yf, cnt, cf):
        col = CONF_CMAP(cnorm(conf))
        lum = 0.299 * col[0] + 0.587 * col[1] + 0.114 * col[2]  # readable text on the disc
        ax.text(
            x, y, human(c), ha="center", va="center", fontsize=7.5, color="white" if lum < 0.55 else "#222222", zorder=3
        )

    ax.set_xlim(-0.6, ncols - 0.4)
    ax.set_ylim(-0.6, nrows - 0.4)
    ax.invert_yaxis()
    ax.set_aspect("equal")  # keep the discs circular
    ax.set_xticks(range(ncols))
    ax.set_yticks(range(nrows))
    ax.set_xticklabels([a.replace("_", " ") for a in ats], rotation=35, ha="right", fontsize=9.5)
    ax.set_yticklabels([k.replace("_", " ") for k in kls], fontsize=9.5)
    ax.set_xlabel("Agent type (AT)", fontsize=10)
    ax.set_ylabel("Knowledge level (KL)", fontsize=10)
    ax.grid(False)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # standard full-height colorbar for disc color (KLAT confidence). Disc size
    # encodes edge count, and each disc is labeled with its count, so there is no
    # separate size legend.
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("KL x AT weight (confidence proxy)", fontsize=9)
    _save(fig, out, dpi)


def _save(fig, out: Path, dpi: int):
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    png = out.with_suffix(".png")
    if png != out:
        fig.savefig(png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}" + (f" and {png}" if png != out else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metagraph", type=Path, help="Path to metagraph JSON")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("klat.pdf"),
        help="Output base path; _bars/_heatmap suffixes are added",
    )
    ap.add_argument(
        "--mode", choices=["bars", "heatmap", "both"], default="both", help="Which view(s) to render (default: both)"
    )
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()
    # write outputs next to this script (not the cwd) unless an absolute path is given
    if not args.output.is_absolute():
        args.output = Path(__file__).resolve().parent / args.output

    meta = json.loads(args.metagraph.read_text())
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42  # embed TrueType, not Type 3 (journals reject Type 3)
    plt.rcParams["ps.fonttype"] = 42

    stem, suf = args.output.with_suffix(""), args.output.suffix or ".pdf"
    if args.mode in ("bars", "both"):
        make_bars(meta, Path(f"{stem}_bars{suf}"), args.dpi)
    if args.mode in ("heatmap", "both"):
        for scale in ("linear", "sqrt", "log"):
            make_heatmap(meta, Path(f"{stem}_heatmap_{scale}{suf}"), args.dpi, scale)


if __name__ == "__main__":
    main()
