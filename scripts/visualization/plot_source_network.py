#!/usr/bin/env python3
"""
Results-section figure for the KRAKEN paper: a force-directed network of PRIMARY
KNOWLEDGE SOURCE adjacency -- i.e. how often edges from two different primary
sources meet at a shared node (a "handoff", the basis for cross-source traversal).

The script scans the edges JSONL directly (and the nodes JSONL, to type each
shared entity by its Biolink family). It does NOT read the metagraph.

  node          = a primary knowledge source (a neutral rounded label "pill")
  node size     ~ number of edges from that source (the label font scales with it)
  edge          = two sources with a shared-entity handoff (gently curved)
  edge color    = Biolink family of the shared entities that dominate the handoff
                  (same palette as the chord figure). Multi-category shared nodes
                  split their credit across families, so the dominant color is the
                  aggregate winner (see scan()).
  edge width    ~ handoff strength, scaled linearly against the 90th-percentile
                  edge (robust to a single huge outlier; outliers clip at max).
                  Strength is set by --measure: 'distinct' shared entities
                  (hub-proof, default) or 'volume' of adjacent edge pairs
                  (hub-dominated).

The (slow) scan can be cached with --cache: the first run scans the JSONL and
writes the cache; later runs load it and skip the scan entirely (so you can
re-style/re-threshold instantly, and the JSONL paths are no longer needed).

Usage:
    # first run: scan the graph and cache the result (cache lands beside this script)
    python plot_source_network.py EDGES.jsonl --nodes NODES.jsonl --cache scan4.pkl -o net.pdf

    # later runs: reuse the cache by bare name (no JSONL needed), tweak backbone/threshold
    python plot_source_network.py --cache scan4.pkl --top-k 3 --min-handoff 100 -o net.pdf
"""
import argparse
import pickle
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.patches import FancyArrowPatch
from networkx.drawing.nx_pydot import graphviz_layout

from kraken.utils.constants import (
    EDGE_OBJECT,
    EDGE_PRIMARY_KS,
    EDGE_SUBJECT,
    NODE_CATEGORIES,
    NODE_ID,
)
from kraken.utils.kg_io import stream_edges_from_jsonl, stream_nodes_from_jsonl

# --- family palette (kept in sync with plot_metagraph_chord.py) ---------------
FAMILY_COLORS = {
    "purple": "#7E5EB0",
    "blue": "#3E7CB1",
    "green": "#3E9B78",
    "yellow": "#D9A441",
    "pink": "#D24559",
    "mauve": "#9E7580",
    "brown": "#9E6B42",
    "slate": "#6E8E88",
    "gray": "#97A0AD",
}
FAMILY_LABELS = {
    "purple": "Genes & variants",
    "blue": "Proteins",
    "green": "Chemicals & metabolites",
    "yellow": "Pathways & processes",
    "pink": "Diseases & phenotypes",
    "mauve": "Anatomy & cells",
    "brown": "Organisms & taxa",
    "slate": "Procedures & activities",
    "gray": "Other",
}
FAMILY_ORDER = list(FAMILY_LABELS)

SCRIPT_DIR = Path(__file__).resolve().parent  # caches + outputs default here

# node "pills": the label font (hence pill size) scales with a source's edge count,
# so visual weight tracks importance rather than name length.
LABEL_MIN_FONT, LABEL_MAX_FONT = 7.0, 13.0
CHAR_W = 0.60  # avg glyph advance as a fraction of font size (bold sans, approx)
LINE_H = 1.15  # line height as a fraction of font size
PILL_PAD_FRAC = 0.35  # pill padding per side, as a fraction of font size (~ bbox pad)
LABEL_GAP = 7.0  # extra points kept between neighboring pills
EDGE_ALPHA = 0.45  # edge transparency
# a lone thin swatch reads lighter than the figure's thick/overlapping edges, so the
# legend gets a slightly higher alpha + thickness to *look* like the drawn edges
LEGEND_ALPHA = 0.62

# display-only friendly labels for a few infores CURIEs whose identifiers are opaque
# (the underlying infores IDs and all data are unchanged; this only affects pill text).
LABEL_OVERRIDES = {
    "multiomics-multiomics": "multiomics-kg",
    "multiomics-microbiome": "microbiome-kg",
    "multiomics-drugapprovals": "drug-approvals-kg",
    "multiomics-clinicaltrials": "clinical-trials-kg",
}


def category_family(category: str) -> str:
    """Biolink category -> family (mirrors plot_metagraph_chord.category_family)."""
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
    if has(
        "organismtaxon",
        "organismalentity",
        "cellularorganism",
        "cellline",
        "individualorganism",
        "populationofindividualorganisms",
        "cohort",
        "human",
    ):
        return "brown"
    if has("biologicalentity"):
        return "blue"
    if has("food", "chemical", "molecularmixture"):
        return "green"
    return "gray"


def short(ks: str) -> str:
    return ks.replace("infores:", "")


def wrap_label(name, width=14):
    """Break a long source name onto two lines at the most central hyphen."""
    if len(name) <= width or "-" not in name:
        return name
    idxs = [i for i, ch in enumerate(name) if ch == "-"]
    i = min(idxs, key=lambda x: abs(x - len(name) / 2))
    return name[: i + 1] + "\n" + name[i + 1 :]


def topk_backbone(G, k):
    """Keep each node's k strongest handoffs (their union) -- a readable backbone
    for a near-complete weighted graph, and a layout that shows real structure."""
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    for n in G.nodes():
        for m, d in sorted(G[n].items(), key=lambda kv: -kv[1]["weight"])[:k]:
            H.add_edge(n, m, weight=d["weight"])
    return H


def separate_labels(pos, nodes, hw, hh, gap=7.0, iters=2000):
    """Push apart node LABEL boxes so text never collides -- rectangle separation
    in point-space (the axes are sized so 1 data unit == 1 pt). Each overlapping
    pair is pushed along its axis of least penetration, which keeps the layout
    compact while guaranteeing non-overlapping pills."""
    P = {n: np.array(pos[n], dtype=float) for n in nodes}
    for _ in range(iters):
        moved = False
        for a, b in combinations(nodes, 2):
            dx = P[b][0] - P[a][0]
            dy = P[b][1] - P[a][1]
            ox = (hw[a] + hw[b] + gap) - abs(dx)  # x-overlap of the two boxes
            oy = (hh[a] + hh[b] + gap) - abs(dy)  # y-overlap
            if ox > 0 and oy > 0:  # boxes intersect
                if ox <= oy:  # least penetration along x
                    s = ox / 2.0 * (1.0 if dx >= 0 else -1.0)
                    P[a][0] -= s
                    P[b][0] += s
                else:  # least penetration along y
                    s = oy / 2.0 * (1.0 if dy >= 0 else -1.0)
                    P[a][1] -= s
                    P[b][1] += s
                moved = True
        if not moved:
            break
    return {n: (float(P[n][0]), float(P[n][1])) for n in nodes}


def scan(edges_file: Path, nodes_file: Path):
    """Scan nodes (for each entity's Biolink family) then edges, and return
    (src_edges, pair_fam): the per-source edge counts, and for each source pair
    the handoff weight split by the Biolink family of the *shared* entity."""
    # A node can carry several Biolink categories (and gene/protein entities are
    # deliberately conflated, so a node may be both). Rather than force one family
    # via an arbitrary tie-break, each node's credit is split EQUALLY across the
    # DISTINCT families it belongs to; the aggregate over all shared nodes then
    # decides each edge's dominant color. Single-family nodes (the vast majority)
    # are stored as a bare string to keep 14.7M entries light.
    print("Scanning nodes for Biolink families...")
    node_family = {}  # node_id -> family str (single) OR {family: weight} (multi)
    for node in stream_nodes_from_jsonl(nodes_file):
        fams = {category_family(c) for c in node[NODE_CATEGORIES]} or {"gray"}
        node_family[node[NODE_ID]] = next(iter(fams)) if len(fams) == 1 else {f: 1.0 / len(fams) for f in fams}

    print("Scanning edges for source adjacency...")
    node_sources = defaultdict(Counter)  # node -> Counter(source)
    src_edges = Counter()  # source -> #edges
    for e in stream_edges_from_jsonl(edges_file):
        ks = short(e.get(EDGE_PRIMARY_KS, "unknown"))
        s, o = e[EDGE_SUBJECT], e[EDGE_OBJECT]
        src_edges[ks] += 1
        node_sources[s][ks] += 1
        node_sources[o][ks] += 1

    # cross-source handoffs split by the shared entity's Biolink family, under two
    # measures per (pair, family):
    #   volume   = sum(count_a * count_b)  -- adjacency volume (hub-dominated)
    #   distinct = # distinct shared nodes -- breadth of overlap (hub-proof)
    pair_vol = defaultdict(Counter)
    pair_dist = defaultdict(Counter)
    for node, srcs in node_sources.items():
        if len(srcs) < 2:
            continue
        fw = node_family.get(node, "gray")
        fam_wts = fw.items() if isinstance(fw, dict) else ((fw, 1.0),)  # (family, weight)
        for (s1, c1), (s2, c2) in combinations(srcs.items(), 2):
            p = tuple(sorted((s1, s2)))
            for fam, wt in fam_wts:
                pair_vol[p][fam] += c1 * c2 * wt
                pair_dist[p][fam] += wt  # fractional; sums to the true distinct-node count
    return src_edges, dict(pair_vol), dict(pair_dist)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "edges", type=Path, nargs="?", default=None, help="Path to edges JSONL (omit if a valid --cache already exists)"
    )
    ap.add_argument(
        "--nodes", type=Path, default=None, help="Path to nodes JSONL (enables Biolink-family node coloring)"
    )
    ap.add_argument(
        "--min-handoff",
        type=int,
        default=1000,
        help="Floor: ignore source pairs below this handoff count (default 1000)",
    )
    ap.add_argument(
        "--top-k", type=int, default=3, help="Draw each source's K strongest handoffs (backbone; default 3)"
    )
    ap.add_argument(
        "--measure",
        choices=["distinct", "volume"],
        default="distinct",
        help="Handoff strength: 'distinct' shared entities (hub-proof, default) "
        "or 'volume' of adjacent edge pairs (hub-dominated)",
    )
    ap.add_argument("-o", "--output", type=Path, default=Path("source_network.pdf"))
    ap.add_argument(
        "--edge-curve", type=float, default=0.1, help="Edge curvature (arc3 rad); 0 = straight lines (default 0.1)"
    )
    ap.add_argument(
        "--edge-scale",
        choices=["p90", "sqrt", "log"],
        default="p90",
        help="Edge-width scaling: 'p90' (linear vs 90th percentile, outliers "
        "clip; default), 'sqrt' (sqrt(w) vs max), or 'log' (log(1+w) vs max). "
        "sqrt/log compress the range progressively more",
    )
    ap.add_argument(
        "--edge-bow",
        choices=["uniform", "radial"],
        default="radial",
        help="Bow direction: 'radial' (default; each edge bows away from the "
        "layout center, opening the core) or 'uniform' (all same handedness)",
    )
    ap.add_argument(
        "--dock-edges",
        action="store_true",
        help="Dock edges at each pill's boundary (not its hidden center) so an "
        "attaching edge visibly terminates at the pill vs. passing behind it",
    )
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument(
        "--cache",
        type=Path,
        default=None,
        help="Pickle cache of the (slow) scan. A bare NAME (no directory) "
        "resolves next to this script (scripts/visualization/); a path "
        "is used as given.",
    )
    args = ap.parse_args()
    if not args.output.is_absolute():
        args.output = SCRIPT_DIR / args.output
    # caches live alongside this script: a bare --cache NAME (no directory component)
    # resolves there, so e.g. `--cache scan4.pkl` works from any working directory.
    if args.cache is not None and args.cache.parent == Path("."):
        args.cache = SCRIPT_DIR / args.cache.name

    if args.cache and args.cache.exists():
        print(f"Loading cached scan from {args.cache}")
        src_edges, pair_vol, pair_dist = pickle.loads(args.cache.read_bytes())
    else:
        if args.edges is None or args.nodes is None:
            raise SystemExit("Provide the edges AND --nodes JSONL, or a --cache that already exists.")
        src_edges, pair_vol, pair_dist = scan(args.edges, args.nodes)
        if args.cache:
            args.cache.write_bytes(pickle.dumps((src_edges, pair_vol, pair_dist)))
            print(f"Cached scan to {args.cache}")

    pair_fam = pair_dist if args.measure == "distinct" else pair_vol
    pair = {p: sum(fc.values()) for p, fc in pair_fam.items()}  # total handoff strength per pair

    gfull = nx.Graph()
    for (a, b), w in pair.items():
        if w >= args.min_handoff:
            gfull.add_edge(a, b, weight=w)
    if gfull.number_of_nodes() == 0:
        raise SystemExit(f"No source pairs reach --min-handoff={args.min_handoff}.")
    G = topk_backbone(gfull, args.top_k)
    print(
        f"Network: {G.number_of_nodes()} sources | {gfull.number_of_edges()} handoffs "
        f">= {args.min_handoff} | {G.number_of_edges()} drawn (top-{args.top_k} per source)"
    )

    def edge_family(u, v):
        fc = pair_fam.get(tuple(sorted((u, v))))
        return fc.most_common(1)[0][0] if fc else "gray"

    nodes = list(G.nodes())
    mx_e = max(src_edges[n] for n in nodes)
    imp = {n: np.log1p(src_edges[n]) / np.log1p(mx_e) for n in nodes}  # 0..1 importance
    fs = {n: LABEL_MIN_FONT + (LABEL_MAX_FONT - LABEL_MIN_FONT) * imp[n] for n in nodes}
    # display label: apply friendly overrides, then hyphenate underscores (display only)
    lbl = {n: wrap_label(LABEL_OVERRIDES.get(short(n), short(n)).replace("_", "-")) for n in nodes}
    # edge-width normalization refs: 90th percentile (for 'p90') and max (for 'log')
    weights = np.array([d["weight"] for *_, d in G.edges(data=True)], dtype=float)
    ref_w = float(np.percentile(weights, 90)) or 1.0
    max_w = float(weights.max()) or 1.0

    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]

    # label-box half-extents in POINTS (approx, proportional bold sans). Everything
    # below works in point-space, and the axes are sized so 1 data unit == 1 pt, so
    # these estimates line up with the drawn pills.
    hw, hh = {}, {}
    for n in nodes:
        lines = lbl[n].split("\n")
        pad = PILL_PAD_FRAC * fs[n]
        hw[n] = max(len(s) for s in lines) * CHAR_W * fs[n] / 2 + pad
        hh[n] = len(lines) * LINE_H * fs[n] / 2 + pad

    # neato layout for topology, rescaled into point-space, then label-box separation
    raw = graphviz_layout(G, prog="neato")
    xs = [raw[n][0] for n in nodes]
    ys = [raw[n][1] for n in nodes]
    span = max(max(xs) - min(xs), max(ys) - min(ys)) or 1.0
    spread = 85.0 * np.sqrt(len(nodes))  # initial layout span in points
    x0, y0 = min(xs), min(ys)
    pos = {n: ((raw[n][0] - x0) / span * spread, (raw[n][1] - y0) / span * spread) for n in nodes}
    pos = separate_labels(pos, nodes, hw, hh, gap=LABEL_GAP)

    # bounds include the label boxes; size the figure so 1 data unit == 1 point
    xmin = min(pos[n][0] - hw[n] for n in nodes) - 28
    xmax = max(pos[n][0] + hw[n] for n in nodes) + 28
    ymin = min(pos[n][1] - hh[n] for n in nodes) - 28
    ymax = max(pos[n][1] + hh[n] for n in nodes) + 28
    fig, ax = plt.subplots(figsize=((xmax - xmin) / 72.0, (ymax - ymin) / 72.0))
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_rasterization_zorder(2)  # flatten curved edges (zorder 1) to raster; text stays vector

    # edge widths: 'p90' = linear vs 90th pct (wider bulk, outliers clip at max width);
    # 'log' = log(1+w) vs the true max (compresses range, lifts small edges, no clip)
    def width(w):
        if args.edge_scale == "log":
            return 1.5 + 6.0 * np.log1p(w) / np.log1p(max_w)
        if args.edge_scale == "sqrt":
            return 1.5 + 6.0 * (w / max_w) ** 0.5
        return 1.5 + 6.0 * min(w / ref_w, 1.0)

    edge_lw = {(u, v): width(d["weight"]) for u, v, d in G.edges(data=True)}

    def dock(n, ux, uy):
        # distance from a pill's center to its box boundary along unit dir (ux, uy),
        # so a docked edge stops just at the pill edge (+small margin) instead of center
        tx = hw[n] / abs(ux) if ux else 1e9
        ty = hh[n] / abs(uy) if uy else 1e9
        return min(tx, ty) + 1.5

    # edges colored by the dominant shared-entity family; curvature from --edge-curve.
    # 'radial' bow flips each edge's arc sign so its midpoint bulges AWAY from the
    # layout centroid -- edges arc around the dense core instead of cutting through it.
    cx = float(np.mean([pos[n][0] for n in nodes]))
    cy = float(np.mean([pos[n][1] for n in nodes]))

    def bow_rad(u, v):
        if args.edge_curve == 0 or args.edge_bow == "uniform":
            return args.edge_curve
        (x0, y0), (x1, y1) = pos[u], pos[v]
        # cross((v-u), (centroid - midpoint)) > 0 => centroid is left of u->v, where
        # a positive rad bows; negate there so the bow goes the other way (outward).
        cross = (x1 - x0) * (cy - (y0 + y1) / 2) - (y1 - y0) * (cx - (x0 + x1) / 2)
        return -args.edge_curve if cross > 0 else args.edge_curve

    for u, v in G.edges():
        sA = sB = 0.0
        if args.dock_edges:
            dx, dy = pos[v][0] - pos[u][0], pos[v][1] - pos[u][1]
            L = (dx * dx + dy * dy) ** 0.5 or 1.0
            sA, sB = dock(u, dx / L, dy / L), dock(v, dx / L, dy / L)
        ax.add_patch(
            FancyArrowPatch(
                pos[u],
                pos[v],
                connectionstyle=f"arc3,rad={bow_rad(u, v)}",
                arrowstyle="-",
                color=FAMILY_COLORS[edge_family(u, v)],
                alpha=EDGE_ALPHA,
                lw=edge_lw[(u, v)],
                capstyle="round",
                shrinkA=sA,
                shrinkB=sB,
                zorder=1,
            )
        )

    # node "pills": a flat, borderless rounded label chip lifted off the edges by a
    # soft shadow. Font size tracks edge count; weight adds only a light hierarchy
    # (regular satellites -> medium hubs, never bold -- the pill does the separating).
    shadow = pe.withSimplePatchShadow(offset=(1.1, -1.3), shadow_rgbFace=(0, 0, 0), alpha=0.16)
    for n in nodes:
        t = ax.text(
            pos[n][0],
            pos[n][1],
            lbl[n],
            ha="center",
            va="center",
            zorder=3,
            fontsize=fs[n],
            fontweight=int(400 + 140 * imp[n]),
            color="#24272C",
            linespacing=0.95,
            bbox=dict(boxstyle="round,pad=0.34", facecolor="#F4F5F7", edgecolor="none"),
        )
        t.get_bbox_patch().set_path_effects([shadow])

    edge_fams = {edge_family(u, v) for u, v in G.edges()}
    present = [f for f in FAMILY_ORDER if f in edge_fams]
    handles = [
        plt.Line2D(
            [], [], color=FAMILY_COLORS[f], lw=6.0, alpha=LEGEND_ALPHA, label=FAMILY_LABELS[f], solid_capstyle="round"
        )
        for f in present
    ]
    ax.legend(
        handles=handles,
        title="Shared node type family",
        loc="upper left",
        frameon=False,
        fontsize=10.5,
        title_fontsize=11.5,
        borderaxespad=1.2,
    )

    fig.savefig(args.output, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    png = args.output.with_suffix(".png")
    if png != args.output:
        fig.savefig(png, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    print(f"Wrote {args.output} and {png}")


if __name__ == "__main__":
    main()
