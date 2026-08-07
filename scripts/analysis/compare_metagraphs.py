"""
Knowledge Graph Meta Graph Comparison Tool

Compares meta graphs for specified sources
to analyze overlapping vs distinct sources, category/predicate coverage,
and help inform graph selection decisions.
"""

import json
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path
from typing import Any


@dataclass
class KGMetaGraph:
    """Parsed representation of a KG meta graph."""

    name: str
    total_nodes: int = 0
    total_edges: int = 0
    node_categories: dict[str, int] = field(default_factory=dict)
    edge_predicates: dict[str, int] = field(default_factory=dict)
    knowledge_sources: dict[str, int] = field(default_factory=dict)
    meta_doubles: dict[str, int] = field(default_factory=dict)
    meta_triples: dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, name: str, data: dict) -> "KGMetaGraph":
        summary = data.get("summary", {})
        return cls(
            name=name,
            total_nodes=summary.get("total_nodes", sum(data.get("node_categories", {}).values())),
            total_edges=summary.get("total_edges", sum(data.get("edge_predicates", {}).values())),
            node_categories=data.get("node_categories", {}),
            edge_predicates=data.get("edge_predicates", {}),
            # metagraphs now split provenance by role; compare on primary sources
            # (fall back to the legacy combined key for older metagraph files)
            knowledge_sources=data.get("primary_knowledge_sources", data.get("knowledge_sources", {})),
            meta_doubles=data.get("meta_doubles", {}),
            meta_triples=data.get("meta_triples", {}),
        )


@dataclass
class KGComparison:
    """Results of comparing multiple KG meta graphs."""

    graphs: dict[str, KGMetaGraph]

    # Source analysis
    @property
    def all_sources(self) -> set[str]:
        return set().union(*(set(g.knowledge_sources.keys()) for g in self.graphs.values()))

    @property
    def shared_sources(self) -> set[str]:
        source_sets = [set(g.knowledge_sources.keys()) for g in self.graphs.values()]
        return set.intersection(*source_sets) if source_sets else set()

    def sources_unique_to(self, kg_name: str) -> set[str]:
        if kg_name not in self.graphs:
            return set()
        this_sources = set(self.graphs[kg_name].knowledge_sources.keys())
        other_sources = set().union(
            *(set(g.knowledge_sources.keys()) for name, g in self.graphs.items() if name != kg_name)
        )
        return this_sources - other_sources

    def sources_shared_by_pair(self, kg1: str, kg2: str) -> set[str]:
        """Get sources shared by exactly these two KGs (and possibly others)."""
        if kg1 not in self.graphs or kg2 not in self.graphs:
            return set()
        sources1 = set(self.graphs[kg1].knowledge_sources.keys())
        sources2 = set(self.graphs[kg2].knowledge_sources.keys())
        return sources1 & sources2

    def sources_exclusive_to_pair(self, kg1: str, kg2: str) -> set[str]:
        """Get sources shared by exactly these two KGs and no others."""
        if kg1 not in self.graphs or kg2 not in self.graphs:
            return set()
        shared = self.sources_shared_by_pair(kg1, kg2)
        other_names = [n for n in self.graphs.keys() if n not in (kg1, kg2)]
        other_sources = (
            set().union(*(set(self.graphs[n].knowledge_sources.keys()) for n in other_names)) if other_names else set()
        )
        return shared - other_sources

    # Category analysis
    @property
    def all_categories(self) -> set[str]:
        return set().union(*(set(g.node_categories.keys()) for g in self.graphs.values()))

    @property
    def shared_categories(self) -> set[str]:
        cat_sets = [set(g.node_categories.keys()) for g in self.graphs.values()]
        return set.intersection(*cat_sets) if cat_sets else set()

    def categories_unique_to(self, kg_name: str) -> set[str]:
        if kg_name not in self.graphs:
            return set()
        this_cats = set(self.graphs[kg_name].node_categories.keys())
        other_cats = set().union(*(set(g.node_categories.keys()) for name, g in self.graphs.items() if name != kg_name))
        return this_cats - other_cats

    # Predicate analysis
    @property
    def all_predicates(self) -> set[str]:
        return set().union(*(set(g.edge_predicates.keys()) for g in self.graphs.values()))

    @property
    def shared_predicates(self) -> set[str]:
        pred_sets = [set(g.edge_predicates.keys()) for g in self.graphs.values()]
        return set.intersection(*pred_sets) if pred_sets else set()

    def predicates_unique_to(self, kg_name: str) -> set[str]:
        if kg_name not in self.graphs:
            return set()
        this_preds = set(self.graphs[kg_name].edge_predicates.keys())
        other_preds = set().union(
            *(set(g.edge_predicates.keys()) for name, g in self.graphs.items() if name != kg_name)
        )
        return this_preds - other_preds

    # Meta triple analysis
    @property
    def all_meta_triples(self) -> set[tuple[str, str, str]]:
        """Returns set of (subject_cat, predicate, object_cat) tuples."""
        triples = set()
        for g in self.graphs.values():
            for meta_triple_str, count in g.meta_triples.items():
                subj, pred, obj = meta_triple_str.split("__")
                triples.add((subj, pred, obj))
        return triples

    def meta_triples_for(self, kg_name: str) -> set[tuple[str, str, str]]:
        if kg_name not in self.graphs:
            return set()
        triples = set()
        g = self.graphs[kg_name]
        for meta_triple_str, count in g.meta_triples.items():
            subj, pred, obj = meta_triple_str.split("__")
            triples.add((subj, pred, obj))
        return triples

    def meta_triples_unique_to(self, kg_name: str) -> set[tuple[str, str, str]]:
        this_triples = self.meta_triples_for(kg_name)
        other_triples = set().union(*(self.meta_triples_for(name) for name in self.graphs.keys() if name != kg_name))
        return this_triples - other_triples


def compare_kgs(*kg_data: tuple[str, dict]) -> KGComparison:
    """Compare multiple KG meta graphs."""
    graphs = {name: KGMetaGraph.from_dict(name, data) for name, data in kg_data}
    return KGComparison(graphs=graphs)


def print_summary_stats(comparison: KGComparison) -> None:
    """Print basic summary statistics."""
    print("=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)

    headers = ["Metric"] + list(comparison.graphs.keys())
    rows = [
        ["Total Nodes"] + [f"{g.total_nodes:,}" for g in comparison.graphs.values()],
        ["Total Edges"] + [f"{g.total_edges:,}" for g in comparison.graphs.values()],
        ["Node Categories"] + [str(len(g.node_categories)) for g in comparison.graphs.values()],
        ["Edge Predicates"] + [str(len(g.edge_predicates)) for g in comparison.graphs.values()],
        ["Primary Sources"] + [str(len(g.knowledge_sources)) for g in comparison.graphs.values()],
        ["Meta Triples"] + [str(len(comparison.meta_triples_for(name))) for name in comparison.graphs.keys()],
    ]

    col_widths = [max(len(str(row[i])) for row in [headers] + rows) for i in range(len(headers))]

    print(" | ".join(h.ljust(w) for h, w in zip(headers, col_widths)))
    print("-+-".join("-" * w for w in col_widths))
    for row in rows:
        print(" | ".join(str(c).ljust(w) for c, w in zip(row, col_widths)))
    print()


def print_source_analysis(comparison: KGComparison) -> None:
    """Print detailed source overlap analysis."""
    print("=" * 70)
    print("SOURCE ANALYSIS")
    print("=" * 70)

    print(f"\nTotal unique sources across all KGs: {len(comparison.all_sources)}")
    print(f"Sources shared by ALL KGs: {len(comparison.shared_sources)}")

    if comparison.shared_sources:
        print("\nShared sources (all KGs):")
        for src in sorted(comparison.shared_sources):
            counts = [f"{name}: {g.knowledge_sources.get(src, 0):,}" for name, g in comparison.graphs.items()]
            print(f"  {src}: {', '.join(counts)}")

    # Pairwise analysis (only when 3+ KGs)
    kg_names = list(comparison.graphs.keys())
    if len(kg_names) > 2:
        print("\n" + "-" * 70)
        print("PAIRWISE SOURCE OVERLAP")
        print("-" * 70)

        for kg1, kg2 in combinations(kg_names, 2):
            shared = comparison.sources_shared_by_pair(kg1, kg2)
            exclusive = comparison.sources_exclusive_to_pair(kg1, kg2)

            print(f"\n{kg1} ∩ {kg2}:")
            print(f"  Total shared: {len(shared)}")
            print(f"  Exclusive to this pair (not in others): {len(exclusive)}")

            if exclusive:
                print("  Exclusive sources:")
                g1, g2 = comparison.graphs[kg1], comparison.graphs[kg2]
                for src in sorted(
                    exclusive,
                    key=lambda s: g1.knowledge_sources.get(s, 0) + g2.knowledge_sources.get(s, 0),
                    reverse=True,
                )[:10]:
                    c1 = g1.knowledge_sources.get(src, 0)
                    c2 = g2.knowledge_sources.get(src, 0)
                    print(f"    {src}: {kg1}={c1:,}, {kg2}={c2:,}")
                if len(exclusive) > 10:
                    print(f"    ... and {len(exclusive) - 10} more")

    print("\nUnique sources per KG:")
    for name in comparison.graphs.keys():
        unique = comparison.sources_unique_to(name)
        print(f"\n  {name} ({len(unique)} unique):")
        if unique:
            g = comparison.graphs[name]
            for src in sorted(unique, key=lambda s: g.knowledge_sources.get(s, 0), reverse=True)[:15]:
                print(f"    {src}: {g.knowledge_sources[src]:,} edges")
            if len(unique) > 15:
                print(f"    ... and {len(unique) - 15} more")
        else:
            print("    (none)")
    print()


def print_category_analysis(comparison: KGComparison) -> None:
    """Print category coverage analysis."""
    print("=" * 70)
    print("NODE CATEGORY ANALYSIS")
    print("=" * 70)

    print(f"\nTotal unique categories: {len(comparison.all_categories)}")
    print(f"Shared by all KGs: {len(comparison.shared_categories)}")

    print("\nCategory coverage (top 20 by max count):")
    all_cats = comparison.all_categories
    cat_max = {cat: max(g.node_categories.get(cat, 0) for g in comparison.graphs.values()) for cat in all_cats}
    top_cats = sorted(all_cats, key=lambda c: cat_max[c], reverse=True)[:20]

    print(f"\n{'Category':<40} " + " ".join(f"{name:>12}" for name in comparison.graphs.keys()))
    print("-" * (40 + 13 * len(comparison.graphs)))
    for cat in top_cats:
        counts = [comparison.graphs[name].node_categories.get(cat, 0) for name in comparison.graphs.keys()]
        count_strs = [f"{c:>12,}" if c > 0 else f"{'--':>12}" for c in counts]
        print(f"{cat:<40} " + " ".join(count_strs))

    print("\nUnique categories per KG:")
    for name in comparison.graphs.keys():
        unique = comparison.categories_unique_to(name)
        if unique:
            print(f"\n  {name}: {sorted(unique)}")
    print()


def print_predicate_analysis(comparison: KGComparison) -> None:
    """Print predicate coverage analysis."""
    print("=" * 70)
    print("EDGE PREDICATE ANALYSIS")
    print("=" * 70)

    print(f"\nTotal unique predicates: {len(comparison.all_predicates)}")
    print(f"Shared by all KGs: {len(comparison.shared_predicates)}")

    print("\nPredicate coverage (top 20 by max count):")
    all_preds = comparison.all_predicates
    pred_max = {pred: max(g.edge_predicates.get(pred, 0) for g in comparison.graphs.values()) for pred in all_preds}
    top_preds = sorted(all_preds, key=lambda p: pred_max[p], reverse=True)[:20]

    print(f"\n{'Predicate':<45} " + " ".join(f"{name:>12}" for name in comparison.graphs.keys()))
    print("-" * (45 + 13 * len(comparison.graphs)))
    for pred in top_preds:
        counts = [comparison.graphs[name].edge_predicates.get(pred, 0) for name in comparison.graphs.keys()]
        count_strs = [f"{c:>12,}" if c > 0 else f"{'--':>12}" for c in counts]
        print(f"{pred:<45} " + " ".join(count_strs))

    print("\nUnique predicates per KG:")
    for name in comparison.graphs.keys():
        unique = comparison.predicates_unique_to(name)
        if unique:
            print(f"\n  {name}: {sorted(unique)}")
    print()


def print_meta_triple_analysis(comparison: KGComparison) -> None:
    """Print meta triple (schema-level path) analysis."""
    print("=" * 70)
    print("META TRIPLE ANALYSIS (Schema-level paths)")
    print("=" * 70)

    all_triples = comparison.all_meta_triples
    print(f"\nTotal unique meta triples: {len(all_triples)}")

    triple_sets = [comparison.meta_triples_for(name) for name in comparison.graphs.keys()]
    shared = set.intersection(*triple_sets) if triple_sets else set()
    print(f"Shared by all KGs: {len(shared)}")

    print("\nUnique meta triples per KG:")
    for name in comparison.graphs.keys():
        unique = comparison.meta_triples_unique_to(name)
        print(f"\n  {name} ({len(unique)} unique meta triples):")
        if unique:
            for subj, pred, obj in sorted(unique)[:10]:
                short_pred = pred.replace("biolink:", "")
                short_subj = subj.replace("biolink:", "")
                short_obj = obj.replace("biolink:", "")
                print(f"    {short_subj} --[{short_pred}]--> {short_obj}")
            if len(unique) > 10:
                print(f"    ... and {len(unique) - 10} more")
    print()


def recommend_kg_for_task(comparison: KGComparison) -> None:
    """Print recommendations for which KG to use for different tasks."""
    print("=" * 70)
    print("KG SELECTION RECOMMENDATIONS")
    print("=" * 70)

    print("\nBest KG for specific node types (by count):")
    categories_of_interest = [
        "biolink:Gene",
        "biolink:Disease",
        "biolink:ChemicalEntity",
        "biolink:Protein",
        "biolink:Pathway",
        "biolink:PhenotypicFeature",
        "biolink:SequenceVariant",
        "biolink:Drug",
        "biolink:AnatomicalEntity",
    ]

    for cat in categories_of_interest:
        counts = {name: g.node_categories.get(cat, 0) for name, g in comparison.graphs.items()}
        if any(counts.values()):
            best = max(counts, key=counts.get)
            print(f"  {cat.replace('biolink:', '')}: {best} ({counts[best]:,})")

    print("\nBest KG for specific relationship types (by count):")
    predicates_of_interest = [
        "biolink:treats",
        "biolink:interacts_with",
        "biolink:associated_with",
        "biolink:causes",
        "biolink:regulates",
        "biolink:expressed_in",
        "biolink:has_phenotype",
        "biolink:gene_associated_with_condition",
    ]

    for pred in predicates_of_interest:
        counts = {name: g.edge_predicates.get(pred, 0) for name, g in comparison.graphs.items()}
        if any(counts.values()):
            best = max(counts, key=counts.get)
            print(f"  {pred.replace('biolink:', '')}: {best} ({counts[best]:,})")

    print("\nSource diversity:")
    for name, g in comparison.graphs.items():
        unique_sources = len(comparison.sources_unique_to(name))
        total_sources = len(g.knowledge_sources)
        print(f"  {name}: {total_sources} sources ({unique_sources} unique)")
    print()


def generate_comparison_report(comparison: KGComparison) -> str:
    """Generate a structured comparison report as a string."""
    import io
    import sys

    old_stdout = sys.stdout
    sys.stdout = buffer = io.StringIO()

    print_summary_stats(comparison)
    print_source_analysis(comparison)
    print_category_analysis(comparison)
    print_predicate_analysis(comparison)
    print_meta_triple_analysis(comparison)
    recommend_kg_for_task(comparison)

    output = buffer.getvalue()
    sys.stdout = old_stdout
    return output


def export_comparison_data(comparison: KGComparison) -> dict[str, Any]:
    """Export comparison data as a dictionary for further analysis."""
    kg_names = list(comparison.graphs.keys())

    # Build pairwise source data
    pairwise_sources = {}
    if len(kg_names) > 2:
        for kg1, kg2 in combinations(kg_names, 2):
            pair_key = f"{kg1}+{kg2}"
            pairwise_sources[pair_key] = {
                "shared": sorted(comparison.sources_shared_by_pair(kg1, kg2)),
                "exclusive": sorted(comparison.sources_exclusive_to_pair(kg1, kg2)),
            }

    return {
        "summary": {
            name: {
                "total_nodes": g.total_nodes,
                "total_edges": g.total_edges,
                "num_categories": len(g.node_categories),
                "num_predicates": len(g.edge_predicates),
                "num_sources": len(g.knowledge_sources),
                "num_meta_triples": len(comparison.meta_triples_for(name)),
            }
            for name, g in comparison.graphs.items()
        },
        "sources": {
            "all": sorted(comparison.all_sources),
            "shared": sorted(comparison.shared_sources),
            "unique_per_kg": {name: sorted(comparison.sources_unique_to(name)) for name in comparison.graphs.keys()},
            "pairwise": pairwise_sources,
        },
        "categories": {
            "all": sorted(comparison.all_categories),
            "shared": sorted(comparison.shared_categories),
            "unique_per_kg": {name: sorted(comparison.categories_unique_to(name)) for name in comparison.graphs.keys()},
        },
        "predicates": {
            "all": sorted(comparison.all_predicates),
            "shared": sorted(comparison.shared_predicates),
            "unique_per_kg": {name: sorted(comparison.predicates_unique_to(name)) for name in comparison.graphs.keys()},
        },
        "meta_triples": {
            "total_unique": len(comparison.all_meta_triples),
            "unique_per_kg": {name: len(comparison.meta_triples_unique_to(name)) for name in comparison.graphs.keys()},
        },
    }


def export_source_matrix_tsv(comparison: KGComparison, path: str | Path) -> None:
    """Export a TSV matrix of knowledge sources across KGs.

    Columns: source_name, kg1_edges, kg2_edges, ...
    Each row is a knowledge source, with edge counts (0 if absent).
    """
    path = Path(path)
    kg_names = list(comparison.graphs.keys())
    all_sources = sorted(comparison.all_sources)

    with open(path, "w") as f:
        # Header
        f.write("knowledge_source\t" + "\t".join(kg_names) + "\n")

        # Data rows
        for src in all_sources:
            counts = [str(comparison.graphs[name].knowledge_sources.get(src, 0)) for name in kg_names]
            f.write(f"{src}\t" + "\t".join(counts) + "\n")


def load_meta_graph(path: str | Path) -> tuple[str, dict]:
    """Load a meta graph from a JSON file, deriving name from filename or 'source' field."""
    path = Path(path)
    with open(path) as f:
        data = json.load(f)

    name = data.get("source", path.stem).upper()
    return name, data


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Compare knowledge graph meta graphs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s spoke_meta.json kg2_meta.json robokop_meta.json
  %(prog)s --export comparison.json spoke_meta.json kg2_meta.json
  %(prog)s --export-sources-tsv sources.tsv spoke_meta.json kg2_meta.json robokop_meta.json
  %(prog)s --names SPOKE KG2 ROBOKOP -- spoke.json kg2.json robokop.json
        """,
    )
    parser.add_argument("metagraphs", nargs="+", help="Paths to meta graph JSON files")
    parser.add_argument("--names", nargs="+", help="Custom names for the KGs (must match number of input files)")
    parser.add_argument("--export", metavar="FILE", help="Export structured comparison data to JSON file")
    parser.add_argument(
        "--export-sources-tsv",
        metavar="FILE",
        help="Export source coverage matrix as TSV (sources × KGs with edge counts)",
    )
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress printed report (useful with --export)")

    args = parser.parse_args()

    kg_data = []
    for i, path in enumerate(args.metagraphs):
        name, data = load_meta_graph(path)
        if args.names and i < len(args.names):
            name = args.names[i]
        kg_data.append((name, data))

    if len(kg_data) < 2:
        parser.error("Need at least 2 meta graphs to compare")

    comparison = compare_kgs(*kg_data)

    if not args.quiet:
        report = generate_comparison_report(comparison)
        print(report)

    if args.export:
        export_data = export_comparison_data(comparison)
        with open(args.export, "w") as f:
            json.dump(export_data, f, indent=2)
        print(f"Exported comparison data to {args.export}")

    if args.export_sources_tsv:
        export_source_matrix_tsv(comparison, args.export_sources_tsv)
        print(f"Exported source matrix to {args.export_sources_tsv}")


if __name__ == "__main__":
    main()
