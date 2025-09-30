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


def load_and_process_summary_data(file_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
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

    return result_df, latest_df


def create_delta_barchart(baseline_df: pd.DataFrame, comparison_df: pd.DataFrame,
                         baseline_raw_df: pd.DataFrame, comparison_raw_df: pd.DataFrame,
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

    # Create ordered labels (first by entity type, then by dataset)
    ordered_labels = []
    ordered_deltas = []
    ordered_counts = []  # Store count information for each bar
    colors = []
    entity_type_positions = []  # Track where each entity type ends

    current_position = 0
    for entity_type in entity_type_order:
        entity_start = current_position
        for dataset in dataset_order:
            mask = (merged_df['dataset'] == dataset) & (merged_df['entity_type'] == entity_type)
            if mask.any():
                row = merged_df[mask].iloc[0]
                entity_display = entity_type_mapping.get(entity_type, entity_type.title())
                # Special formatting for dataset names
                if dataset.lower() == 'ukbb':
                    dataset_display = 'UKBB'
                elif dataset.lower() == 'israeli10k':
                    dataset_display = 'HPP'
                else:
                    dataset_display = dataset.title()
                label = f"{dataset_display}"  # Just the dataset name
                ordered_labels.append(label)
                ordered_deltas.append(row['delta'])

                # Store count information for display - get actual counts from raw data
                baseline_raw_row = baseline_raw_df[(baseline_raw_df['dataset'] == dataset) &
                                                  (baseline_raw_df['entity_type'] == entity_type)]
                comparison_raw_row = comparison_raw_df[(comparison_raw_df['dataset'] == dataset) &
                                                      (comparison_raw_df['entity_type'] == entity_type)]

                if len(baseline_raw_row) > 0 and len(comparison_raw_row) > 0:
                    # Get actual mapped counts and total items
                    baseline_mapped = int(baseline_raw_row['mapped_to_kg'].iloc[0])
                    comparison_mapped = int(comparison_raw_row['mapped_to_kg'].iloc[0])
                    total_items = int(comparison_raw_row['total_items'].iloc[0])

                    # Calculate the difference
                    count_diff = comparison_mapped - baseline_mapped

                    # Always show count info now
                    count_info = f"({count_diff:+d} / {total_items})"
                else:
                    count_info = ""

                ordered_counts.append(count_info)

                current_position += 1
                # Color coding: positive = green, negative = red, zero = gray
                if row['delta'] > 0:
                    colors.append('#2E8B57')  # Sea green
                elif row['delta'] < 0:
                    colors.append('#DC143C')  # Crimson
                else:
                    colors.append('#808080')  # Gray

        # Record the end position of this entity type (if it had any bars)
        if current_position > entity_start:
            entity_type_positions.append(current_position)

    # Create the bar chart
    plt.figure(figsize=(10, 8))

    bars = plt.bar(range(len(ordered_labels)), ordered_deltas, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5, width=0.6)

    # Customize the plot
    plt.title(f'Coverage Delta: {comparison_name} vs {baseline_name}', fontsize=16, fontweight='bold', pad=50)
    plt.xlabel('Dataset', fontsize=12, fontweight='bold')
    plt.ylabel('Coverage Difference (%)', fontsize=12, fontweight='bold')

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

    # Add vertical lines to separate entity type groups
    for i in range(len(entity_type_positions) - 1):  # Don't draw after the last group
        x_position = entity_type_positions[i] - 0.5
        plt.axvline(x=x_position, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    # Add entity type labels right at the top of the chart area
    current_pos = 0
    for i, entity_type in enumerate(entity_type_order):
        # Find the range for this entity type
        if i < len(entity_type_positions):
            end_pos = entity_type_positions[i]
            if end_pos > current_pos:  # Only add label if there are bars for this entity type
                center_pos = (current_pos + end_pos - 1) / 2
                entity_display = entity_type_mapping.get(entity_type, entity_type.title())
                plt.text(center_pos, plt.ylim()[1] - 1, entity_display,
                        ha='center', va='bottom', fontsize=11, fontweight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgray', alpha=1.0))
                current_pos = end_pos

    # Add delta labels and count labels above bars
    for i, (bar, delta, counts) in enumerate(zip(bars, ordered_deltas, ordered_counts)):
        height = bar.get_height()
        x_center = bar.get_x() + bar.get_width()/2.

        if height >= 0:  # Positive bars - stack labels above
            # Delta percentage at the top
            delta_y = height + 3
            # Count info below the delta (closer spacing)
            count_y = height + 1
            va_delta = 'bottom'
            va_count = 'bottom'
        else:  # Negative bars - stack labels below
            # Delta percentage at the bottom
            delta_y = height - 3
            # Count info above the delta (closer spacing)
            count_y = height - 1
            va_delta = 'top'
            va_count = 'top'

        # Delta label
        plt.text(x_center, delta_y, f'{delta:+.1f}%',
                ha='center', va=va_delta, fontsize=9, fontweight='bold')

        # Count label (always show now if there's text)
        if counts:  # Only add label if counts is not empty
            plt.text(x_center, count_y, counts,
                    ha='center', va=va_count, fontsize=7, fontweight='normal')

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
    baseline_df, baseline_raw_df = load_and_process_summary_data(baseline_file)

    print(f"Loading comparison data from: {comparison_file}")
    comparison_df, comparison_raw_df = load_and_process_summary_data(comparison_file)

    print(f"\nBaseline data shape: {baseline_df.shape}")
    print(f"Comparison data shape: {comparison_df.shape}")

    # Create visualization
    merged_df = create_delta_barchart(baseline_df, comparison_df, baseline_raw_df, comparison_raw_df,
                                     args.baseline_name, args.comparison_name, output_dir)

    # Generate summary statistics
    generate_delta_summary(merged_df, args.baseline_name, args.comparison_name, output_dir)


if __name__ == "__main__":
    main()