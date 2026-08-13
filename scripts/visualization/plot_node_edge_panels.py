#!/usr/bin/env python3
"""
Compose two example images into a single two-panel figure (node example on the
left, edge example on the right) WITHOUT distorting either image's proportions.

Each panel is given the exact aspect ratio of its source image and both panels
share a common height, so images are only ever scaled uniformly -- never
stretched. Optional bold panel labels (A, B) sit above the panels.

Output is written as PDF (raster images embedded, vector labels) plus a PNG,
next to this script unless an absolute -o path is given.

Usage:
    python plot_node_edge_panels.py NODE.png EDGE.png -o node_edge.pdf
    python plot_node_edge_panels.py NODE.png EDGE.png --height 3.0 --gap 0.3
    python plot_node_edge_panels.py NODE.png EDGE.png --labels ""   # no A/B labels
"""
import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent


def load(path: Path):
    """Load an image and return (array, aspect = width/height)."""
    if not path.exists():
        raise SystemExit(f"Image not found: {path}")
    img = mpimg.imread(path)
    h, w = img.shape[0], img.shape[1]
    return img, w / h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("node_image", type=Path, help="Left panel image (e.g. the node example)")
    ap.add_argument("edge_image", type=Path, help="Right panel image (e.g. the edge example)")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("node_edge_panels.pdf"),
        help="Output path; a matching .png is also written (default beside this script)",
    )
    ap.add_argument("--height", type=float, default=3.2, help="Panel height in inches (default 3.2)")
    ap.add_argument("--gap", type=float, default=0.28, help="Gap between panels in inches (default 0.28)")
    ap.add_argument("--labels", default="A,B", help="Comma-separated panel labels (default 'A,B'; pass '' for none)")
    ap.add_argument("--label-size", type=float, default=6.0, help="Panel-label font size (default 15)")
    ap.add_argument("--dpi", type=int, default=2000)
    ap.add_argument(
        "--background",
        default="white",
        help="Canvas color behind/around the panels (default #FAFAFA; makes white "
        "screenshots pop). Pass 'white' for none.",
    )
    ap.add_argument(
        "--border", default="#E6E6E6", help="Thin border color around each panel (default #E6E6E6; 'none' to omit)"
    )
    ap.add_argument(
        "--pad", type=float, default=0.15, help="Background margin around the whole figure, in inches (default 0.15)"
    )
    args = ap.parse_args()
    if not args.output.is_absolute():
        args.output = SCRIPT_DIR / args.output

    imgL, aspL = load(args.node_image)
    imgR, aspR = load(args.edge_image)
    labels = [s for s in (t.strip() for t in args.labels.split(",")) if s] if args.labels.strip() else []

    plt.rcParams["pdf.fonttype"] = 42  # embed TrueType, not Type 3 (journals reject Type 3)
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]

    # both panels share height H; widths follow each image's aspect -> no distortion
    H = args.height
    wL, wR = H * aspL, H * aspR
    total_w = wL + args.gap + wR
    label_band = 0.34 if labels else 0.0  # inches reserved above the panels for A/B
    total_h = H + label_band
    frac_h = H / total_h  # panels occupy the bottom band of the figure

    fig = plt.figure(figsize=(total_w, total_h), facecolor=args.background)
    # each axes box is placed at exactly its image's aspect ratio, so aspect="auto"
    # fills the box without stretching the pixels
    axL = fig.add_axes([0.0, 0.0, wL / total_w, frac_h])
    axR = fig.add_axes([(wL + args.gap) / total_w, 0.0, wR / total_w, frac_h])
    for ax, img in ((axL, imgL), (axR, imgR)):
        ax.imshow(img, aspect="auto", interpolation="antialiased")
        ax.set_xticks([])
        ax.set_yticks([])
        # keep a thin panel border so the (white) screenshots read against the canvas
        for spine in ax.spines.values():
            if args.border.lower() == "none":
                spine.set_visible(False)
            else:
                spine.set_edgecolor(args.border)
                spine.set_linewidth(0.2)

    if labels:
        y = frac_h + (1.0 - frac_h) * 0.18  # just above the panels
        fig.text(0.0, y, labels[0], ha="left", va="bottom", fontsize=args.label_size, fontweight="bold")
        if len(labels) > 1:
            fig.text(
                (wL + args.gap) / total_w,
                y,
                labels[1],
                ha="left",
                va="bottom",
                fontsize=args.label_size,
                fontweight="bold",
            )

    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight", facecolor=args.background, pad_inches=args.pad)
    png = args.output.with_suffix(".png")
    if png != args.output:
        fig.savefig(png, dpi=args.dpi, bbox_inches="tight", facecolor=args.background, pad_inches=args.pad)
    plt.close(fig)
    print(f"Wrote {args.output}" + (f" and {png}" if png != args.output else ""))


if __name__ == "__main__":
    main()
