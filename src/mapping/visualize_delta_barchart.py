#!/usr/bin/env python3
"""
Visualize delta bar chart comparing mapping coverage between two knowledge graphs.
Shows the difference in coverage percentages for each entity type and dataset.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re


def parse_dataset_info(dataset_name: str) -> tuple[str, str, int]:
    """Parse dataset name into dataset, entity type, and version"""
    # Pattern: {dataset}_{entity_type}_v{version}
    match = re.match(r'([^_]+)_([^_]+)_v(\d+)', dataset_name)
    if match:
        dataset = match.group(1)
        entity_type = match.group(2)
        version = int(match.group(3))
        return dataset, entity_type, version
    else:
        raise ValueError(f"Cannot parse dataset name: {dataset_name}")


def load_and_process_summary_data(file_path: Path) -> pd.DataFrame:
    """Load TSV data and process into coverage DataFrame"""

    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')

    # Remove empty rows
    df = df.dropna(how='all')

    # Parse dataset information and add columns
    parsed_data = []
    for _, row in df.iterrows():
        try:
            dataset, entity_type, version = parse_dataset_info(row['dataset'])
            parsed_data.append({
                'original_dataset': row['dataset'],
                'dataset': dataset,
                'entity_type': entity_type,
                'version': version,
                'total_items': row['total_items'],
                'mapped_to_kg': row['mapped_to_kg']
            })
        except ValueError as e:
            print(f"Warning: {e}")
            continue

    parsed_df = pd.DataFrame(parsed_data)

    # Keep only the latest version for each dataset/entity_type pair
    latest_versions = parsed_df.groupby(['dataset', 'entity_type'])['version'].max().reset_index()
    latest_df = parsed_df.merge(latest_versions, on=['dataset', 'entity_type', 'version'])

    # Calculate coverage percentages
    latest_df['coverage_percentage'] = (latest_df['mapped_to_kg'] / latest_df['total_items']) * 100

    # Create a matrix with dataset and entity type as index
    result_df = latest_df.pivot_table(
        index=['dataset', 'entity_type'],
        values='coverage_percentage',
        aggfunc='first'
    ).reset_index()

    return result_df


def create_delta_barchart(baseline_df: pd.DataFrame, comparison_df: pd.DataFrame,
                         baseline_name: str, comparison_name: str, output_dir: Path):
    """Create and save the delta bar chart"""

    # Merge the dataframes to calculate deltas
    merged_df = baseline_df.merge(
        comparison_df,
        on=['dataset', 'entity_type'],
        how='outer',
        suffixes=('_baseline', '_comparison')
    )

    # Fill NaN values with 0 for calculation
    merged_df['coverage_percentage_baseline'] = merged_df['coverage_percentage_baseline'].fillna(0)
    merged_df['coverage_percentage_comparison'] = merged_df['coverage_percentage_comparison'].fillna(0)

    # Calculate delta (comparison - baseline)
    merged_df['delta'] = merged_df['coverage_percentage_comparison'] - merged_df['coverage_percentage_baseline']

    # Create a combined label for dataset_entitytype
    merged_df['dataset_entity'] = merged_df['dataset'] + '_' + merged_df['entity_type']

    # Order datasets: arivale, ukbb, israeli10k
    dataset_order = ['arivale', 'ukbb', 'israeli10k']

    # Map entity type names for display
    entity_type_mapping = {
        'clinicallabs': 'Clinical Labs',
        'lipids': 'Lipids',
        'metabolites': 'Metabolites',
        'proteins': 'Proteins'
    }

    # Order entity types in reverse alphabetical order (same as heatmap)
    entity_type_order = ['proteins', 'metabolites', 'lipids', 'clinicallabs']

    # Create ordered labels
    ordered_labels = []
    ordered_deltas = []
    colors = []

    for dataset in dataset_order:
        for entity_type in entity_type_order:
            mask = (merged_df['dataset'] == dataset) & (merged_df['entity_type'] == entity_type)
            if mask.any():
                row = merged_df[mask].iloc[0]
                entity_display = entity_type_mapping.get(entity_type, entity_type.title())
                label = f"{dataset.title()}\n{entity_display}"
                ordered_labels.append(label)
                ordered_deltas.append(row['delta'])
                # Color coding: positive = green, negative = red, zero = gray
                if row['delta'] > 0:
                    colors.append('#2E8B57')  # Sea green
                elif row['delta'] < 0:
                    colors.append('#DC143C')  # Crimson
                else:
                    colors.append('#808080')  # Gray

    # Create the bar chart
    plt.figure(figsize=(10, 8))

    bars = plt.bar(range(len(ordered_labels)), ordered_deltas, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)

    # Customize the plot
    plt.title(f'Coverage Delta: {comparison_name} vs {baseline_name}', fontsize=16, fontweight='bold', pad=50)
    plt.xlabel('Dataset - Entity Type', fontsize=12, fontweight='bold')
    plt.ylabel('Coverage Difference (percent)', fontsize=12, fontweight='bold')

    # Set x-axis labels
    plt.xticks(range(len(ordered_labels)), ordered_labels, rotation=45, ha='right')

    # Set consistent y-axis scale for comparability
    min_delta = min(ordered_deltas) if ordered_deltas else 0
    max_delta = max(ordered_deltas) if ordered_deltas else 0

    if min_delta >= 0:  # Only positive or zero differences
        plt.ylim(0, 100)
    elif max_delta <= 0:  # Only negative or zero differences
        plt.ylim(-100, 0)
    else:  # Both positive and negative differences
        plt.ylim(-100, 100)

    # Add horizontal line at y=0
    plt.axhline(y=0, color='black', linestyle='-', linewidth=1, alpha=0.5)

    # Add value labels on bars
    for i, (bar, delta) in enumerate(zip(bars, ordered_deltas)):
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + (1 if height >= 0 else -3),
                f'{delta:+.1f}%', ha='center', va='bottom' if height >= 0 else 'top',
                fontsize=9, fontweight='bold')

    # Add legend outside the plot area
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#2E8B57', alpha=0.7, label=f'{comparison_name} gain vs. {baseline_name}'),
        Patch(facecolor='#DC143C', alpha=0.7, label=f'{comparison_name} loss vs. {baseline_name}')
    ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(0.5, 1.04), loc='lower center', ncol=2)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # Make room for rotated labels and footnote

    # Save the plot
    output_file = output_dir / f'delta_coverage_{comparison_name.lower()}_vs_{baseline_name.lower()}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved delta bar chart: {output_file}")

    # Also save as PDF
    output_file_pdf = output_dir / f'delta_coverage_{comparison_name.lower()}_vs_{baseline_name.lower()}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved delta bar chart (PDF): {output_file_pdf}")

    plt.close()

    return merged_df


def generate_delta_summary(merged_df: pd.DataFrame, baseline_name: str, comparison_name: str, output_dir: Path):
    """Generate and save delta summary statistics"""

    print(f"\n=== DELTA COVERAGE STATISTICS: {comparison_name} vs {baseline_name} ===")

    # Overall statistics
    deltas = merged_df['delta'].dropna()
    if not deltas.empty:
        mean_delta = deltas.mean()
        max_delta = deltas.max()
        min_delta = deltas.min()

        print(f"Average delta: {mean_delta:+.1f} percentage points")
        print(f"Maximum delta: {max_delta:+.1f} percentage points")
        print(f"Minimum delta: {min_delta:+.1f} percentage points")

        # Count improvements/degradations
        improvements = (deltas > 0).sum()
        degradations = (deltas < 0).sum()
        equal = (deltas == 0).sum()

        print(f"\nSummary:")
        print(f"  {comparison_name} better: {improvements} cases")
        print(f"  {baseline_name} better: {degradations} cases")
        print(f"  Equal performance: {equal} cases")

        # Dataset-wise statistics
        print(f"\nDataset-wise average delta:")
        for dataset in merged_df['dataset'].unique():
            if pd.notna(dataset):
                dataset_deltas = merged_df[merged_df['dataset'] == dataset]['delta'].dropna()
                if not dataset_deltas.empty:
                    avg_delta = dataset_deltas.mean()
                    print(f"  {dataset}: {avg_delta:+.1f} percentage points")

        # Entity type-wise statistics
        print(f"\nEntity type-wise average delta:")
        for entity_type in merged_df['entity_type'].unique():
            if pd.notna(entity_type):
                entity_deltas = merged_df[merged_df['entity_type'] == entity_type]['delta'].dropna()
                if not entity_deltas.empty:
                    avg_delta = entity_deltas.mean()
                    entity_display = {'clinicallabs': 'Clinical Labs', 'lipids': 'Lipids',
                                    'metabolites': 'Metabolites', 'proteins': 'Proteins'}.get(entity_type, entity_type.title())
                    print(f"  {entity_display}: {avg_delta:+.1f} percentage points")

    # Save summary to file
    summary_file = output_dir / f'delta_coverage_stats_{comparison_name.lower()}_vs_{baseline_name.lower()}.txt'
    with open(summary_file, 'w') as f:
        f.write(f"=== DELTA COVERAGE STATISTICS: {comparison_name} vs {baseline_name} ===\n\n")

        if not deltas.empty:
            f.write(f"Average delta: {mean_delta:+.1f} percentage points\n")
            f.write(f"Maximum delta: {max_delta:+.1f} percentage points\n")
            f.write(f"Minimum delta: {min_delta:+.1f} percentage points\n\n")

            f.write("Summary:\n")
            f.write(f"  {comparison_name} better: {improvements} cases\n")
            f.write(f"  {baseline_name} better: {degradations} cases\n")
            f.write(f"  Equal performance: {equal} cases\n\n")

            f.write("Dataset-wise average delta:\n")
            for dataset in merged_df['dataset'].unique():
                if pd.notna(dataset):
                    dataset_deltas = merged_df[merged_df['dataset'] == dataset]['delta'].dropna()
                    if not dataset_deltas.empty:
                        avg_delta = dataset_deltas.mean()
                        f.write(f"  {dataset}: {avg_delta:+.1f} percentage points\n")

            f.write(f"\nEntity type-wise average delta:\n")
            for entity_type in merged_df['entity_type'].unique():
                if pd.notna(entity_type):
                    entity_deltas = merged_df[merged_df['entity_type'] == entity_type]['delta'].dropna()
                    if not entity_deltas.empty:
                        avg_delta = entity_deltas.mean()
                        entity_display = {'clinicallabs': 'Clinical Labs', 'lipids': 'Lipids',
                                        'metabolites': 'Metabolites', 'proteins': 'Proteins'}.get(entity_type, entity_type.title())
                        f.write(f"  {entity_display}: {avg_delta:+.1f} percentage points\n")

    print(f"Saved delta summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description='Create delta bar chart comparing coverage between two knowledge graphs')
    parser.add_argument('baseline_file', help='Baseline summary TSV file path')
    parser.add_argument('comparison_file', help='Comparison summary TSV file path')
    parser.add_argument('--baseline-name', default='BASELINE', help='Name for baseline KG (default: BASELINE)')
    parser.add_argument('--comparison-name', default='COMPARISON', help='Name for comparison KG (default: COMPARISON)')
    parser.add_argument('--output', default='results', help='Output directory (default: src/mapping)')
    args = parser.parse_args()

    # File paths
    baseline_file = Path(args.baseline_file)
    comparison_file = Path(args.comparison_file)
    output_dir = Path(args.output)

    if not baseline_file.exists():
        print(f"Error: Baseline file {baseline_file} not found!")
        return

    if not comparison_file.exists():
        print(f"Error: Comparison file {comparison_file} not found!")
        return

    # Load and process data
    print(f"Loading baseline data from: {baseline_file}")
    baseline_df = load_and_process_summary_data(baseline_file)

    print(f"Loading comparison data from: {comparison_file}")
    comparison_df = load_and_process_summary_data(comparison_file)

    print(f"\nBaseline data shape: {baseline_df.shape}")
    print(f"Comparison data shape: {comparison_df.shape}")

    # Create visualization
    merged_df = create_delta_barchart(baseline_df, comparison_df, args.baseline_name, args.comparison_name, output_dir)

    # Generate summary statistics
    generate_delta_summary(merged_df, args.baseline_name, args.comparison_name, output_dir)

    print(f"\nUsage tips:")
    print(f"  python src/mapping/visualize_delta_barchart.py baseline.tsv comparison.tsv")
    print(f"  python src/mapping/visualize_delta_barchart.py kraken.tsv kg2.tsv --baseline-name KRAKEN --comparison-name KG2")


if __name__ == "__main__":
    main()