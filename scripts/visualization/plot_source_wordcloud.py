#!/usr/bin/env python3
"""
Graphical-abstract word cloud of KRAKEN's primary knowledge sources, sized by how
many edges each source contributes (log-scaled by default so the long tail of
smaller sources stays legible next to the giants).

A few sources can be excluded by substring match on their infores CURIE via
--exclude (default: pgs-catalog, nih-cde, loinc -- placeholders we plan to add,
so the filter is future-proofed and simply no-ops until they appear).

Colors are drawn from the same Biolink "family" palette as the other figures for
a consistent look (or pass --colormap to use a matplotlib colormap instead).

Outputs a high-resolution PNG and a vector SVG (crisp at any size) next to this
script unless an absolute -o is given.

Usage:
    python plot_source_wordcloud.py METAGRAPH.json -o source_wordcloud.png
    python plot_source_wordcloud.py METAGRAPH.json --scale sqrt --background none
"""
import argparse
import json
import zlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex, to_rgb
from PIL import Image
from wordcloud import WordCloud

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_FONT = "/System/Library/Fonts/HelveticaNeue.ttc"

# same family palette as the network/chord figures (gray dropped -- too faint here)
PALETTE = ["#7E5EB0", "#3E7CB1", "#3E9B78", "#D9A441", "#D24559", "#9E7580", "#9E6B42", "#6E8E88"]


def _mute(hexc, f=0.42):
    """Blend a color toward mid-gray (softer, less saturated)."""
    return to_hex(tuple(v * (1 - f) + 0.5 * f for v in to_rgb(hexc)))


MUTED = [_mute(c) for c in PALETTE]
STEEL = LinearSegmentedColormap.from_list("steel", ["#B7C7D6", "#5E86A6", "#26485F"])  # size gradient

TRANSFORM = {"log": np.log1p, "sqrt": np.sqrt, "linear": lambda x: np.asarray(x, float)}


def _pick(palette, word):
    """Deterministic per-word color from a palette (stable across runs)."""
    return palette[zlib.crc32(word.encode()) % len(palette)]


def make_mask(shape, w, h):
    """Binary WordCloud mask (255 = no words) for a non-box layout region."""
    yy, xx = np.ogrid[:h, :w]
    cx, cy = w / 2.0, h / 2.0
    if shape == "circle":
        r = min(w, h) / 2.0
        inside = ((xx - cx) / r) ** 2 + ((yy - cy) / r) ** 2 <= 1
    elif shape == "blob":  # superellipse: a soft rounded rectangle
        inside = np.abs((xx - cx) / (w / 2.0)) ** 3.2 + np.abs((yy - cy) / (h / 2.0)) ** 3.2 <= 1
    else:  # ellipse fills the full width/height
        inside = ((xx - cx) / (w / 2.0)) ** 2 + ((yy - cy) / (h / 2.0)) ** 2 <= 1
    return np.where(inside, 0, 255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metagraph", type=Path, help="Path to metagraph JSON")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("source_wordcloud.png"),
        help="Output image; a matching .svg is also written (default beside this script)",
    )
    ap.add_argument(
        "--exclude",
        default="pgs-catalog,nih-cde,loinc",
        help="Comma-separated substrings; any source whose infores CURIE contains one "
        "is dropped (default: pgs-catalog,nih-cde,loinc)",
    )
    ap.add_argument(
        "--scale",
        choices=["log", "sqrt", "linear"],
        default="log",
        help="How edge counts map to word size (default log)",
    )
    ap.add_argument(
        "--width", type=int, default=2400, help="Canvas width px (roomy so all words fit at the capped size range)"
    )
    ap.add_argument("--height", type=int, default=1350)
    ap.add_argument(
        "--background", default="transparent", help="Background color, or 'none'/'transparent' for a transparent PNG"
    )
    ap.add_argument(
        "--color-mode",
        choices=["family", "muted", "mono", "size", "cmap"],
        default="size",
        help="Coloring: 'family' palette, 'muted' (desaturated family), 'mono' (one "
        "color), 'size' (single-hue gradient, darker=more edges), or 'cmap' "
        "(random from --colormap). Default family",
    )
    ap.add_argument("--mono-color", default="#33373E", help="Word color for --color-mode mono")
    ap.add_argument(
        "--colormap", default=None, help="Matplotlib colormap for --color-mode size/cmap (defaults to a steel ramp)"
    )
    ap.add_argument(
        "--shape",
        choices=["box", "ellipse", "circle", "blob"],
        default="blob",
        help="Built-in layout region shape (default box)",
    )
    ap.add_argument(
        "--mask", type=Path, default=None, help="Custom mask image (black shape on white); overrides --shape"
    )
    ap.add_argument("--font", default=DEFAULT_FONT, help="Path to a .ttf/.ttc font")
    ap.add_argument(
        "--prefer-horizontal",
        type=float,
        default=0.9,
        help="Fraction of words laid out horizontally (1.0 = all horizontal)",
    )
    ap.add_argument(
        "--size-ratio",
        type=float,
        default=2.5,
        help="Cap on the biggest word's size vs the smallest: font size is made "
        "proportional to the (log) edge count, remapped so max/min = this ratio "
        "(default 2.5; lower = more uniform)",
    )
    ap.add_argument(
        "--max-font-size",
        type=float,
        default=80.0,
        help="Size (px) of the biggest word; the smallest is this / --size-ratio",
    )
    ap.add_argument("--max-words", type=int, default=250)
    ap.add_argument("--png-scale", type=float, default=2.0, help="Upscale factor for the raster PNG")
    ap.add_argument("--seed", type=int, default=42, help="Layout random seed (reproducible)")
    args = ap.parse_args()
    if not args.output.is_absolute():
        args.output = SCRIPT_DIR / args.output

    pks = json.loads(args.metagraph.read_text())["primary_knowledge_sources"]
    excludes = [s.strip() for s in args.exclude.split(",") if s.strip()]

    weights = {}
    dropped = []
    for curie, count in pks.items():
        if any(x in curie for x in excludes):
            dropped.append(curie)
            continue
        name = curie[len("infores:") :] if curie.startswith("infores:") else curie
        weights[name.replace("_", "-")] = float(TRANSFORM[args.scale]([count])[0])

    # remap the (log) weights into [1/size_ratio, 1] so that -- with relative_scaling=1
    # below -- font size is proportional to the weight and the max/min ratio is exactly
    # size_ratio. This is the honest scaling knob; WordCloud's own sizing is not.
    lo, hi = min(weights.values()), max(weights.values())
    span = (hi - lo) or 1.0
    floor = 1.0 / args.size_ratio
    freqs = {n: floor + (1.0 - floor) * (v - lo) / span for n, v in weights.items()}

    if dropped:
        print(f"excluded {len(dropped)}: {', '.join(dropped)}")
    print(f"{len(freqs)} sources in the cloud (of {len(pks)} total)")

    transparent = args.background.lower() in ("none", "transparent", "")

    mask = None
    if args.mask is not None:
        mask = np.array(Image.open(args.mask).convert("L"))
        mask = np.where(mask > 200, 255, 0).astype(np.uint8)  # near-white -> masked out
    elif args.shape != "box":
        mask = make_mask(args.shape, args.width, args.height)

    # static per-word color function by mode (size is applied after layout, below)
    if args.color_mode == "muted":
        color_func = lambda word, *a, **k: _pick(MUTED, word)  # noqa: E731
    elif args.color_mode == "mono":
        color_func = lambda word, *a, **k: args.mono_color  # noqa: E731
    elif args.color_mode == "cmap":
        color_func = None  # let WordCloud sample --colormap
    else:  # family (also the placeholder for 'size', recolored afterwards)
        color_func = lambda word, *a, **k: _pick(PALETTE, word)  # noqa: E731

    wc = WordCloud(
        width=args.width,
        height=args.height,
        mask=mask,
        background_color=None if transparent else args.background,
        mode="RGBA" if transparent else "RGB",
        font_path=args.font,
        prefer_horizontal=args.prefer_horizontal,
        relative_scaling=1.0,  # font size ∝ our remapped weight
        # clamp BOTH ends so WordCloud can't shrink small words to fit -- it drops them
        # instead, keeping the true max/min ratio at size_ratio (not 14x)
        max_font_size=args.max_font_size,
        min_font_size=max(4, round(args.max_font_size / args.size_ratio)),
        max_words=args.max_words,
        colormap=args.colormap if args.color_mode == "cmap" else None,
        color_func=color_func,
        random_state=args.seed,
        margin=4,
        collocations=False,
        scale=args.png_scale,
    )
    wc.generate_from_frequencies(freqs)

    sizes = [item[1] for item in wc.layout_]
    print(
        f"placed {len(sizes)}/{len(freqs)} words; font {min(sizes)}-{max(sizes)}px "
        f"(ratio {max(sizes) / max(min(sizes), 1):.2f})"
    )

    if args.color_mode == "size":  # color each word by its font size (darker = more edges)
        norm = Normalize(min(sizes), max(sizes))
        cmap = plt.get_cmap(args.colormap) if args.colormap else STEEL
        wc.recolor(color_func=lambda word, font_size, *a, **k: to_hex(cmap(norm(font_size))))

    wc.to_file(str(args.output))
    svg = args.output.with_suffix(".svg")
    try:
        svg.write_text(wc.to_svg(embed_font=True))
    except Exception:  # .ttc collections can't be embedded; reference the font by name instead
        svg.write_text(wc.to_svg(embed_font=False))
    print(f"Wrote {args.output} and {svg}")


if __name__ == "__main__":
    main()
