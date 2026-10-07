#!/usr/bin/env python3
"""
Build the KRAKEN supplementary-tables workbook (.xlsx) directly from a metagraph
JSON. A README/cover sheet (graph version, summary stats, notes) plus four simple
two-column, count-sorted tables:

  S1. Primary knowledge sources  (infores CURIE, edge count)    <- primary_knowledge_sources
  S2. Vocabulary prefixes        (prefix, node count)           <- node_prefixes
  S3. Node types                 (Biolink category, node count) <- node_categories
  S4. Edge types                 (Biolink predicate, edge count)<- edge_predicates

Optionally appends two further sections (each only when its input is given):

  S5. Example node and edge records  <- --example-node / --example-edge (one JSON record each)
  S6. Figure queries                 <- --figure-queries (JSON: {"note": ..., "queries": [{figure, description, url}]})

Usage:
    python make_supplementary_tables.py METAGRAPH.json -o supplementary_tables.xlsx
    python make_supplementary_tables.py METAGRAPH.json -o supplementary_data.xlsx \
        --example-node supplementary_example_node.json --example-edge supplementary_example_edge.json \
        --figure-queries supplementary_figure_queries.json
"""
import argparse
import json
from datetime import date
from pathlib import Path

from openpyxl import Workbook
from openpyxl.styles import Alignment, Font

SCRIPT_DIR = Path(__file__).resolve().parent

AGGREGATOR_NOTE = (
    "Note: robokop-kg, rtx-kg2, translator-kg-open and sri-node-normalizer are aggregator or identifier-mapping "
    "resources, but also appear here because some edges are attributed directly to them rather than to an "
    "upstream source (for example close_match cross-references and, for robokop-kg, member_of edges)."
)

MULTI_NOTE = (
    "Note: these counts sum to MORE than the total number of KRAKEN nodes, because a "
    "single node can carry multiple {}. See the README sheet."
)

TABLES = [
    dict(
        title="S1. Primary knowledge sources",
        key="primary_knowledge_sources",
        kh="Primary knowledge source (infores CURIE)",
        vh="Edge count",
        desc="Every primary knowledge source and the number of edges attributed to it.",
        note=AGGREGATOR_NOTE,
    ),
    dict(
        title="S2. Vocabulary prefixes",
        key="node_prefixes",
        kh="Vocabulary prefix",
        vh="Node count",
        desc="Every identifier prefix and the number of node identifiers that use it.",
        note=MULTI_NOTE.format("equivalent identifiers (and thus multiple prefixes)"),
    ),
    dict(
        title="S3. Node types",
        key="node_categories",
        kh="Node type (Biolink category)",
        vh="Node count",
        desc="Every node type and the number of nodes assigned to it.",
        note=MULTI_NOTE.format("Biolink categories"),
    ),
    dict(
        title="S4. Edge types",
        key="edge_predicates",
        kh="Edge type (Biolink predicate)",
        vh="Edge count",
        desc="Every edge type and the number of edges of that type.",
        note=None,
    ),
]

EXAMPLES_SECTION = dict(
    title="S5. Example node and edge records",
    desc="One node record and one edge record as they appear in KRAKEN's NDJSON files.",
)
QUERIES_SECTION = dict(
    title="S6. Figure queries",
    desc="Web interface queries used to generate the manuscript's figures.",
)
EXAMPLES_INTRO = (
    "Each line of KRAKEN's nodes and edges NDJSON files is one JSON record. The records below are the node and "
    "edge shown in Figure 1, reproduced from the files and indented here for readability (in the files each "
    "record occupies a single line)."
)

README_NOTES = [
    "In S2 (prefixes) and S3 (node types), the counts sum to MORE than the total number of nodes: "
    "a single node can carry multiple equivalent identifiers (hence multiple vocabulary prefixes) "
    "and multiple Biolink categories.",
    "In S1 (primary knowledge sources) and S4 (edge types), the counts sum to the total number of "
    "edges: each edge has exactly one primary knowledge source and one predicate.",
    "S1 includes four aggregator or identifier-mapping resources (robokop-kg, rtx-kg2, translator-kg-open and "
    "sri-node-normalizer): some edges are attributed directly to them rather than to an upstream source, for "
    "example close_match cross-references and, for robokop-kg, member_of edges.",
]

# (metagraph summary key, README label) -- only those present are shown
SUMMARY_ROWS = [
    ("total_nodes", "Total nodes"),
    ("total_edges", "Total edges"),
    ("unique_node_categories", "Unique node types (Biolink categories)"),
    ("unique_node_prefixes", "Unique vocabulary prefixes"),
    ("unique_edge_predicates", "Unique edge types (Biolink predicates)"),
    ("unique_primary_knowledge_sources", "Unique primary knowledge sources"),
    ("unique_aggregator_knowledge_sources", "Unique aggregator knowledge sources"),
    ("unique_supporting_data_sources", "Unique supporting data sources"),
]


def _count_cell(cell):
    cell.number_format = "#,##0"
    cell.alignment = Alignment(horizontal="right")


def write_table(ws, t, data: dict):
    row = 1
    if t["note"]:  # short caveat above the header
        ws.cell(row, 1, t["note"])
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=2)
        c = ws.cell(row, 1)
        c.font = Font(italic=True, color="595959")
        c.alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[row].height = 28 if len(t["note"]) < 160 else 56
        row += 1

    ws.cell(row, 1, t["kh"]).font = Font(bold=True)
    ws.cell(row, 2, t["vh"]).font = Font(bold=True)
    ws.freeze_panes = ws.cell(row + 1, 1).coordinate

    max_key = len(t["kh"])
    for name, count in sorted(data.items(), key=lambda kv: -kv[1]):
        row += 1
        ws.cell(row, 1, name)
        _count_cell(ws.cell(row, 2, count))
        max_key = max(max_key, len(str(name)))

    ws.column_dimensions["A"].width = min(max_key + 2, 60)
    ws.column_dimensions["B"].width = max(len(t["vh"]) + 2, 12)


def build_readme(wb, meta, source_name, extras=()):
    ws = wb.create_sheet("README", 0)
    ws.sheet_view.showGridLines = False
    ws.column_dimensions["A"].width = 46
    ws.column_dimensions["B"].width = 40
    r = 1

    def title(text, size):
        nonlocal r
        ws.cell(r, 1, text).font = Font(bold=True, size=size)
        r += 1

    def kv(label, value, count=False):
        nonlocal r
        ws.cell(r, 1, label)
        cell = ws.cell(r, 2, value)
        if count:
            _count_cell(cell)
        r += 1

    def note(text):
        nonlocal r
        ws.cell(r, 1, "• " + text)
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=2)
        ws.cell(r, 1).alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[r].height = 42
        r += 1

    title("KRAKEN knowledge graph: supplementary data", 14)
    r += 1
    kv("Graph name", meta.get("graph", "—"))
    kv("Graph version", meta.get("version", "—"))
    kv("Biolink Model version", meta.get("biolink_version", "—"))
    kv("Source metagraph", source_name)
    kv("Generated", date.today().isoformat())
    r += 1

    title("Summary statistics", 12)
    summ = meta.get("summary", {})
    for key, label in SUMMARY_ROWS:
        if key in summ:
            kv(label, summ[key], count=True)
    r += 1

    title("Contents", 12)
    for t in [*TABLES, *extras]:
        ws.cell(r, 1, t["title"]).font = Font(bold=True)
        ws.cell(r, 2, t["desc"]).alignment = Alignment(wrap_text=True, vertical="top")
        r += 1
    r += 1

    title("Notes", 12)
    for n in README_NOTES:
        note(n)


def _json_lines(record: dict) -> list[str]:
    return json.dumps(record, indent=2, ensure_ascii=False).splitlines()


def write_examples_sheet(wb, node: dict | None, edge: dict | None):
    ws = wb.create_sheet(title="S5. Example records")
    ws.column_dimensions["A"].width = 120
    ws.cell(1, 1, EXAMPLES_INTRO).alignment = Alignment(wrap_text=True, vertical="top")
    ws.row_dimensions[1].height = 44
    r = 3
    for label, record in (("Example node record", node), ("Example edge record", edge)):
        if record is None:
            continue
        ws.cell(r, 1, label).font = Font(bold=True)
        r += 1
        for ln in _json_lines(record):
            ws.cell(r, 1, ln).font = Font(name="Menlo", size=9)
            r += 1
        r += 1


def write_queries_sheet(wb, fq: dict):
    ws = wb.create_sheet(title=QUERIES_SECTION["title"])
    for col, width in (("A", 12), ("B", 80), ("C", 120)):
        ws.column_dimensions[col].width = width
    r = 1
    if fq.get("note"):
        ws.cell(r, 1, fq["note"]).alignment = Alignment(wrap_text=True, vertical="top")
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=3)
        ws.row_dimensions[r].height = 44
        r += 1
    for col, header in enumerate(("Figure", "Description", "Query link"), start=1):
        ws.cell(r, col, header).font = Font(bold=True)
    ws.freeze_panes = ws.cell(r + 1, 1).coordinate
    for q in fq["queries"]:
        r += 1
        ws.cell(r, 1, q["figure"])
        ws.cell(r, 2, q["description"]).alignment = Alignment(wrap_text=True, vertical="top")
        link = ws.cell(r, 3, q["url"])
        if len(q["url"]) <= 2000:  # Excel's hyperlink length limit is ~2079 characters
            link.hyperlink = q["url"]
            link.font = Font(color="0563C1", underline="single")


class _TextPages:
    """Flowing text for the PDF: writes lines top-down and starts a new page when one fills up."""

    def __init__(self, pdf, plt, page, title):
        self.pdf, self.plt, self.page, self.title = pdf, plt, page, title
        self.n = 0
        self.fig = None
        self._new_page()

    def _new_page(self):
        self._flush()
        self.n += 1
        self.fig = self.plt.figure(figsize=self.page)
        self.fig.text(0.08, 0.955, self.title + ("  (continued)" if self.n > 1 else ""), fontsize=13, fontweight="bold")
        self.y = 0.918

    def _flush(self):
        if self.fig is not None:
            self.fig.text(0.5, 0.03, f"{self.title} · page {self.n}", ha="center", fontsize=7, color="#888888")
            self.pdf.savefig(self.fig)
            self.plt.close(self.fig)
            self.fig = None

    def need(self, height):
        """Start a new page unless `height` (figure fraction) still fits on this one."""
        if self.y - height < 0.06:
            self._new_page()

    def line(self, text, size=8.5, mono=False, bold=False, color="black", url=None, x=0.08):
        dy = size * 1.42 / 72 / self.page[1]
        self.need(dy)
        self.fig.text(
            x,
            self.y,
            text,
            fontsize=size,
            family="monospace" if mono else "sans-serif",
            fontweight="bold" if bold else "normal",
            color=color,
            url=url,
            va="top",
        )
        self.y -= dy

    def gap(self, frac=0.012):
        self.y -= frac

    def close(self):
        self._flush()


def _wrap_json_line(ln: str, width: int) -> list[str]:
    """Wrap one pretty-printed JSON line, keeping continuation lines indented under it."""
    import textwrap

    if len(ln) <= width:
        return [ln]
    indent = len(ln) - len(ln.lstrip(" "))
    return textwrap.wrap(ln, width, subsequent_indent=" " * (indent + 4), break_long_words=True) or [ln]


def _pdf_examples(pdf, plt, page, node, edge):
    import textwrap

    tp = _TextPages(pdf, plt, page, EXAMPLES_SECTION["title"])
    for wl in textwrap.wrap(EXAMPLES_INTRO, 118):
        tp.line(wl, size=9)
    for label, record in (("Example node record", node), ("Example edge record", edge)):
        if record is None:
            continue
        wrapped = [wl for ln in _json_lines(record) for wl in _wrap_json_line(ln, 112)]
        line_h = 7.2 * 1.42 / 72 / page[1]
        tp.gap()
        # keep a record on one page when it fits, so it never ends with a stray closing brace overleaf
        body_h = len(wrapped) * line_h
        tp.need(0.035 + body_h if body_h <= 0.8 else 0.08)
        tp.line(label, size=10.5, bold=True)
        tp.gap(0.004)
        for wl in wrapped:
            tp.line(wl, size=7.2, mono=True)
    tp.close()


def _pdf_queries(pdf, plt, page, fq):
    import textwrap

    tp = _TextPages(pdf, plt, page, QUERIES_SECTION["title"])
    for wl in textwrap.wrap(fq.get("note", ""), 118):
        tp.line(wl, size=9)
    for q in fq["queries"]:
        desc = textwrap.wrap(q["description"], 118)
        url_lines = textwrap.wrap(q["url"], 122, break_long_words=True, break_on_hyphens=False)
        tp.gap()
        tp.need(0.02 + 0.014 * (len(desc) + min(len(url_lines), 4)))  # keep heading, description, link start together
        tp.line(q["figure"], size=10.5, bold=True)
        for wl in desc:
            tp.line(wl, size=9)
        tp.gap(0.003)
        # every wrapped line carries the full link, so clicking anywhere on it opens the query
        for wl in url_lines:
            tp.line(wl, size=6.8, mono=True, color="#0563C1", url=q["url"])
    tp.close()


def make_pdf(pdf_path, meta, source_name, node=None, edge=None, figure_queries=None):
    """Render the same content as a paginated, print-ready PDF (title page + tables)."""
    import textwrap

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    plt.rcParams["pdf.fonttype"] = 42  # embed TrueType, not Type 3
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"]
    plt.rcParams["font.monospace"] = ["Menlo", "DejaVu Sans Mono", "Courier New"]
    PAGE = (8.5, 11.0)  # US Letter portrait
    extras = ([EXAMPLES_SECTION] if (node or edge) else []) + ([QUERIES_SECTION] if figure_queries else [])
    ROWS_PER_PAGE = 44

    with PdfPages(pdf_path) as pdf:
        # --- title / README page ---
        fig = plt.figure(figsize=PAGE)
        y = [0.95]

        def line(label, val=None, size=10, bold=False, wrap=None):
            if wrap:
                for wl in textwrap.wrap(label, wrap):
                    fig.text(0.08, y[0], wl, fontsize=size)
                    y[0] -= 0.020
                y[0] -= 0.006
                return
            fig.text(0.08, y[0], label, fontsize=size, fontweight="bold" if bold else "normal")
            if val is not None:
                fig.text(0.62, y[0], val, fontsize=size)
            y[0] -= 0.030

        fig.text(0.08, y[0], "KRAKEN knowledge graph: supplementary data", fontsize=15, fontweight="bold")
        y[0] -= 0.055
        line("Graph name", str(meta.get("graph", "—")))
        line("Graph version", str(meta.get("version", "—")))
        line("Biolink Model version", str(meta.get("biolink_version", "—")))
        line("Source metagraph", source_name)
        line("Generated", date.today().isoformat())
        y[0] -= 0.02
        line("Summary statistics", size=12, bold=True)
        summ = meta.get("summary", {})
        for key, label in SUMMARY_ROWS:
            if key in summ:
                line(label, f"{summ[key]:,}")
        y[0] -= 0.02
        line("Contents", size=12, bold=True)
        for t in [*TABLES, *extras]:
            fig.text(0.08, y[0], t["title"], fontsize=9.5, fontweight="bold")
            fig.text(0.40, y[0], t["desc"], fontsize=9)
            y[0] -= 0.030
        y[0] -= 0.02
        line("Notes", size=12, bold=True)
        for n in README_NOTES:
            line("• " + n, size=9, wrap=105)
        pdf.savefig(fig)
        plt.close(fig)

        # --- one (paginated) section per table ---
        for t in TABLES:
            items = sorted(meta[t["key"]].items(), key=lambda kv: -kv[1])
            chunks = [items[i : i + ROWS_PER_PAGE] for i in range(0, len(items), ROWS_PER_PAGE)] or [[]]
            for pi, chunk in enumerate(chunks):
                fig = plt.figure(figsize=PAGE)
                title = t["title"] + ("  (continued)" if pi else "")
                fig.text(0.08, 0.955, title, fontsize=13, fontweight="bold")
                yh = 0.915
                fig.text(0.08, yh, t["kh"], fontsize=8.5, fontweight="bold")
                fig.text(0.93, yh, t["vh"], fontsize=8.5, fontweight="bold", ha="right")
                fig.add_artist(
                    plt.Line2D(
                        [0.08, 0.93], [yh - 0.008, yh - 0.008], color="#999999", lw=0.6, transform=fig.transFigure
                    )
                )
                dy = 0.855 / ROWS_PER_PAGE
                for ri, (name, count) in enumerate(chunk):
                    yy = yh - 0.022 - ri * dy
                    fig.text(0.08, yy, str(name), fontsize=7.5)
                    fig.text(0.93, yy, f"{count:,}", fontsize=7.5, ha="right")
                fig.text(
                    0.5,
                    0.03,
                    f"{t['title']} · page {pi + 1} of {len(chunks)}",
                    ha="center",
                    fontsize=7,
                    color="#888888",
                )
                pdf.savefig(fig)
                plt.close(fig)

        if node or edge:
            _pdf_examples(pdf, plt, PAGE, node, edge)
        if figure_queries:
            _pdf_queries(pdf, plt, PAGE, figure_queries)

    print(f"Wrote {pdf_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("metagraph", type=Path, help="Path to metagraph JSON")
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("supplementary_tables.xlsx"),
        help="Output .xlsx path (default beside this script)",
    )
    ap.add_argument(
        "--pdf",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Also write a paginated PDF alongside the .xlsx (default: yes; --no-pdf to skip)",
    )
    ap.add_argument("--example-node", type=Path, default=None, help="JSON file holding one example node record (S5)")
    ap.add_argument("--example-edge", type=Path, default=None, help="JSON file holding one example edge record (S5)")
    ap.add_argument(
        "--figure-queries",
        type=Path,
        default=None,
        help="JSON file of the web interface queries behind the figures (S6): "
        '{"note": str, "queries": [{"figure", "description", "url"}]}',
    )
    args = ap.parse_args()
    if not args.output.is_absolute():
        args.output = SCRIPT_DIR / args.output

    meta = json.loads(args.metagraph.read_text())
    node = json.loads(args.example_node.read_text()) if args.example_node else None
    edge = json.loads(args.example_edge.read_text()) if args.example_edge else None
    figure_queries = json.loads(args.figure_queries.read_text()) if args.figure_queries else None
    extras = ([EXAMPLES_SECTION] if (node or edge) else []) + ([QUERIES_SECTION] if figure_queries else [])

    wb = Workbook()
    wb.remove(wb.active)  # drop the default empty sheet

    build_readme(wb, meta, args.metagraph.name, extras)
    for t in TABLES:
        if t["key"] not in meta:
            raise SystemExit(f"metagraph is missing '{t['key']}' (needed for {t['title']})")
        ws = wb.create_sheet(title=t["title"])
        write_table(ws, t, meta[t["key"]])
        print(f"{t['title']}: {len(meta[t['key']]):,} rows")

    if node or edge:
        write_examples_sheet(wb, node, edge)
    if figure_queries:
        write_queries_sheet(wb, figure_queries)

    wb.save(args.output)
    print(f"Wrote {args.output}")

    if args.pdf:
        make_pdf(args.output.with_suffix(".pdf"), meta, args.metagraph.name, node, edge, figure_queries)


if __name__ == "__main__":
    main()
