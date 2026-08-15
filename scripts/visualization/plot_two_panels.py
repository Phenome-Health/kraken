#!/usr/bin/env python3
"""
Compose two images into a two-panel figure WITHOUT distorting either image's
proportions. Vertical (panel A on top, B below) or horizontal (A left, B right).

Panels share one dimension (width when vertical, height when horizontal) and the
other follows each image's aspect ratio, so images are only ever scaled
uniformly. Optional bold panel labels (A, B) sit at each panel's top-left.

Outputs a PDF (raster images embedded, vector labels) plus a PNG, next to this
script unless an absolute -o is given.

Usage:
    python plot_two_panels.py TOP.png BOTTOM.png -o fig.pdf            # vertical (default)
    python plot_two_panels.py LEFT.png RIGHT.png --orientation horizontal
    python plot_two_panels.py A.png B.png --labels ""                  # no A/B labels
"""
import argparse
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent


def load(path: Path):
    """Load an image, returning (array, aspect = width / height)."""
    if not path.exists():
        raise SystemExit(f"Image not found: {path}")
    img = mpimg.imread(path)
    return img, img.shape[1] / img.shape[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("image_a", type=Path, help="Panel A image (top if vertical, left if horizontal)")
    ap.add_argument("image_b", type=Path, help="Panel B image (bottom if vertical, right if horizontal)")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("two_panels.pdf"),
        help="Output path; a matching .png is also written (default beside this script)",
    )
    ap.add_argument(
        "--orientation",
        choices=["vertical", "horizontal"],
        default="vertical",
        help="Stack panels vertically (A over B, default) or horizontally (A | B)",
    )
    ap.add_argument(
        "--size",
        type=float,
        default=6.5,
        help="Shared dimension in inches: panel width if vertical, height if horizontal",
    )
    ap.add_argument("--gap", type=float, default=0.2, help="Gap between panels in inches (default 0.2)")
    ap.add_argument("--labels", default="A,B", help="Comma-separated panel labels (default 'A,B'; pass '' for none)")
    ap.add_argument("--label-size", type=float, default=14.0, help="Panel-label font size (default 14)")
    ap.add_argument("--background", default="white", help="Figure/panel background color (default white)")
    ap.add_argument(
        "--border", default="none", help="Thin border color around each panel ('none' to omit; default none)"
    )
    ap.add_argument("--pad", type=float, default=0.1, help="Margin around the figure, inches (default 0.1)")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()
    if not args.output.is_absolute():
        args.output = SCRIPT_DIR / args.output

    panels = [load(args.image_a), load(args.image_b)]
    labels = [s for s in (t.strip() for t in args.labels.split(",")) if s] if args.labels.strip() else []

    plt.rcParams["pdf.fonttype"] = 42  # embed TrueType, not Type 3 (journals reject Type 3)
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]

    band = 0.30 if labels else 0.0  # inches reserved above each panel for its label
    n = len(panels)

    # per-panel (width, height) in inches; the shared dimension is args.size
    if args.orientation == "vertical":
        dims = [(args.size, args.size / asp) for _, asp in panels]
        total_w = args.size
        total_h = sum(h for _, h in dims) + args.gap * (n - 1) + band * n
    else:  # horizontal
        dims = [(args.size * asp, args.size) for _, asp in panels]
        total_h = args.size + band
        total_w = sum(w for w, _ in dims) + args.gap * (n - 1)

    fig = plt.figure(figsize=(total_w, total_h), facecolor=args.background)

    if args.orientation == "vertical":
        y = 0.0  # inches from the top
        for (img, _), (w, h), i in zip(panels, dims, range(n)):
            panel_top = y + band
            ax = fig.add_axes([0.0, 1 - (panel_top + h) / total_h, w / total_w, h / total_h])
            _draw(ax, img, args)
            if labels:
                fig.text(
                    0.004,
                    1 - panel_top / total_h + 0.006,
                    labels[i],
                    ha="left",
                    va="bottom",
                    fontsize=args.label_size,
                    fontweight="bold",
                )
            y = panel_top + h + args.gap
    else:  # horizontal
        x = 0.0  # inches from the left
        frac_h = args.size / total_h
        for (img, _), (w, h), i in zip(panels, dims, range(n)):
            left = x / total_w
            ax = fig.add_axes([left, 0.0, w / total_w, frac_h])
            _draw(ax, img, args)
            if labels:
                fig.text(
                    left + 0.002,
                    frac_h + (1 - frac_h) * 0.18,
                    labels[i],
                    ha="left",
                    va="bottom",
                    fontsize=args.label_size,
                    fontweight="bold",
                )
            x += w + args.gap

    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight", facecolor=args.background, pad_inches=args.pad)
    png = args.output.with_suffix(".png")
    if png != args.output:
        fig.savefig(png, dpi=args.dpi, bbox_inches="tight", facecolor=args.background, pad_inches=args.pad)
    plt.close(fig)
    print(f"Wrote {args.output}" + (f" and {png}" if png != args.output else ""))


def _draw(ax, img, args):
    """Show an image filling its axes box (box already matches image aspect)."""
    ax.imshow(img, aspect="auto", interpolation="antialiased")
    ax.set_facecolor(args.background)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        if args.border.lower() == "none":
            spine.set_visible(False)
        else:
            spine.set_edgecolor(args.border)
            spine.set_linewidth(0.8)


if __name__ == "__main__":
    main()
