#!/usr/bin/env python3
"""
Figure 1 for the KRAKEN paper: a chord diagram of category-to-category
connectivity, derived from a metagraph JSON file's `meta_doubles`.

Categories are colored by Biolink "family" (adapted from the kestrel-ui
getCategoryColorGroup logic, with organismal entities merged into the taxon
family) and ordered around the circle by family so each family forms a
contiguous block. Connectivity is symmetric (subject/object
directions summed); use the heatmap script for the directed view.

Three versions are written for each run, differing only in how arc/ribbon
widths encode edge counts:
  *_linear : widths proportional to true edge counts (honest proportions).
  *_sqrt   : widths square-root-scaled -- a soft compression that lifts small
             categories while keeping the big ones dominant.
  *_log    : widths log-scaled [log(1+edges)] -- strong compression so even
             low-volume categories are clearly visible.
The sqrt and log versions are disclosed in-figure, because their widths are
NOT proportional to raw counts.

Usage:
    python plot_metagraph_chord.py METAGRAPH.json [-o OUT.pdf] [--top N]
        [--min-ribbon K] [--dpi 300]
"""
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from mpl_chord_diagram import chord_diagram

# --- family palette: deeper than the UI pastels, harmonized with the inventory
# figure (Figure 2). blue matches Figure 2's steel blue so the figures share a spine.
FAMILY_COLORS = {
    "purple": "#7E5EB0",  # genes, variants, transcripts, nucleic acids
    "blue": "#3E7CB1",  # proteins, complexes, biological entity
    "green": "#3E9B78",  # drugs, small molecules, chemicals, food
    "yellow": "#D9A441",  # pathways & processes
    "pink": "#CD5A7C",  # disease, clinical, phenotype
    "mauve": "#A9748A",  # anatomy, cellular component
    "brown": "#9E6B42",  # organism taxa + organismal entities (cell lines, etc.)
    "slate": "#6E8E88",  # procedure, activity, device
    "gray": "#97A0AD",  # everything else
}
# order around the circle (related families adjacent)
FAMILY_ORDER = ["purple", "blue", "green", "yellow", "pink", "mauve", "brown", "slate", "gray"]
FAMILY_LABELS = {
    "purple": "Genes & variants",
    "blue": "Proteins",
    "green": "Chemicals & drugs",
    "yellow": "Pathways & processes",
    "pink": "Diseases & phenotypes",
    "mauve": "Anatomy & cells",
    "brown": "Organisms & taxa",
    "slate": "Procedures & activities",
    "gray": "Other",
}


def category_family(category: str) -> str:
    """Port of kestrel-ui getCategoryColorGroup (same rules, same order)."""
    c = category.replace("biolink:", "").lower()

    def has(*xs):
        return any(x in c for x in xs)

    if has(
        "gene",
        "sequencevariant",
        "snv",
        "transcript",
        "haplotype",
        "microrna",
        "noncodingrnaproduct",
        "genomicentity",
        "nucleicacidentity",
        "rnaproduct",
    ):
        return "purple"
    if has("disease"):
        return "pink"
    if has("clinicalattribute", "clinicalcourse", "clinicalonset", "clinicalmeasurement"):
        return "pink"
    if has("phenotyp", "behavioralfeature", "clinicalfinding"):
        return "pink"
    if has("drug", "smallmolecule", "molecularentity"):
        return "green"
    if has("pathway", "biologicalprocess", "molecularactivity", "pathologicalprocess", "physiologicalprocess"):
        return "yellow"
    if has("macromolecularcomplex", "protein", "polypeptide"):
        return "blue"
    if has("procedure", "activity", "behavior", "device"):
        return "slate"
    if has(
        "anatomicalentity", "cellularcomponent", "cell", "grossanatomicalstructure", "pathologicalanatomicalstructure"
    ):
        return "mauve"
    # organism taxa AND organismal entities (cell lines, individual organisms,
    # populations) are merged into one brown family for the figure -- a deliberate
    # deviation from the UI's getCategoryColorGroup, which keeps them separate.
    if has(
        "organismtaxon",
        "organismalentity",
        "cellularorganism",
        "cellline",
        "individualorganism",
        "populationofindividualorganisms",
    ):
        return "brown"
    if has("biologicalentity"):
        return "blue"
    if has("food", "chemical", "molecularmixture"):
        return "green"
    return "gray"


def short(c: str) -> str:
    return c.replace("biolink:", "")


# width transform per scale. Note: sqrt/log widths are NOT proportional to raw
# edge counts -- disclose that in the figure caption (kept off-figure by request).
SCALES = {
    "linear": lambda M: M,
    "sqrt": np.sqrt,
    "log": np.log1p,
}


def draw(cats, M, colors, families, scale, out, dpi):
    mat = SCALES[scale](M)
    fig, ax = plt.subplots(figsize=(10.5, 10.5), subplot_kw=dict(aspect="equal"))
    chord_diagram(
        mat,
        [short(c) for c in cats],
        ax=ax,
        colors=colors,
        gap=0.03,
        use_gradient=True,
        sort=None,
        chordwidth=0.7,
        fontsize=9,
        rotate_names=True,
    )
    ax.set_xlim(-1.4, 1.4)
    ax.set_ylim(-1.4, 1.4)

    # rasterize the ribbons/arcs (fast, small PDF) but keep text/labels vector.
    # done before the legend so its swatches stay crisp vector.
    for artist in list(ax.patches) + list(ax.collections):
        artist.set_rasterized(True)

    present = [f for f in FAMILY_ORDER if f in set(families)]
    handles = [Patch(facecolor=FAMILY_COLORS[f], label=FAMILY_LABELS[f]) for f in present]
    leg = ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(-0.08, 1.02),
        frameon=False,
        fontsize=8,
        title="Node type family",
        title_fontsize=9,
    )
    leg._legend_box.sep = 10  # extra vertical space under the legend title

    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    png = out.with_suffix(".png")
    if png != out:
        fig.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out}" + (f" and {png}" if png != out else ""))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metagraph", type=Path, help="Path to metagraph JSON")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("chord.pdf"),
        help="Output base path; _linear/_sqrt/_log suffixes are added",
    )
    ap.add_argument("--top", type=int, default=0, help="Keep only the top N categories by connectivity (0 = all)")
    ap.add_argument(
        "--min-ribbon", type=int, default=100, help="Hide ribbons below this raw edge count (declutter; default 100)"
    )
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()
    # write outputs next to this script (not the cwd) unless an absolute path is given
    if not args.output.is_absolute():
        args.output = Path(__file__).resolve().parent / args.output
    print(args)

    meta = json.loads(args.metagraph.read_text())
    doubles = meta["meta_doubles"]

    # rank categories by total connectivity, then order by family
    involve = {}
    for s, objs in doubles.items():
        for o, c in objs.items():
            involve[s] = involve.get(s, 0) + c
            involve[o] = involve.get(o, 0) + c
    cats = sorted(involve, key=lambda c: -involve[c])
    if args.top and args.top > 0:
        cats = cats[: args.top]
    cats.sort(key=lambda c: (FAMILY_ORDER.index(category_family(c)), -involve[c]))

    idx = {c: i for i, c in enumerate(cats)}
    M = np.zeros((len(cats), len(cats)))
    for s, objs in doubles.items():
        if s not in idx:
            continue
        for o, c in objs.items():
            if o in idx:
                M[idx[s], idx[o]] += c
    M = M + M.T
    np.fill_diagonal(M, 0)
    if args.min_ribbon > 0:
        M[M < args.min_ribbon] = 0

    # drop categories left with no off-diagonal connectivity (avoids empty arcs)
    keep = M.sum(axis=1) > 0
    if not keep.all():
        dropped = [cats[i] for i in range(len(cats)) if not keep[i]]
        print(
            f"dropping {len(dropped)} unconnected categories (after applying meta-edge count cutoff): "
            f"{', '.join(short(c) for c in dropped)}"
        )
        cats = [c for i, c in enumerate(cats) if keep[i]]
        M = M[np.ix_(keep, keep)]

    families = [category_family(c) for c in cats]
    colors = [FAMILY_COLORS[f] for f in families]

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    stem, suf = args.output.with_suffix(""), args.output.suffix or ".pdf"
    for scale in SCALES:
        draw(cats, M, colors, families, scale, Path(f"{stem}_{scale}{suf}"), args.dpi)


if __name__ == "__main__":
    main()
