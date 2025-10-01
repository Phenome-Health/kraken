#!/usr/bin/env python3
"""
Visualize mapping summary data as a grid of waterfall charts.
Each cell contains a waterfall chart showing the progression from total items to mapped items.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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
    """Load TSV data and process into structured DataFrame"""

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
                'has_valid_ids': row['has_valid_ids'],
                'has_valid_ids_provided': row['has_valid_ids_provided'],
                'has_only_provided_ids': row['has_only_provided_ids'],
                'has_only_assigned_ids': row['has_only_assigned_ids'],
                'has_both_provided_and_assigned_ids': row['has_both_provided_and_assigned_ids'],
                'mapped_to_kg': row['mapped_to_kg'],
                'one_to_one_mappings': row['one_to_one_mappings'],
                'multi_mappings': row['multi_mappings']
            })
        except ValueError as e:
            print(f"Warning: {e}")
            continue

    parsed_df = pd.DataFrame(parsed_data)

    # Keep only the latest version for each dataset/entity_type pair
    latest_versions = parsed_df.groupby(['dataset', 'entity_type'])['version'].max().reset_index()
    latest_df = parsed_df.merge(latest_versions, on=['dataset', 'entity_type', 'version'])

    print(f"Original entries: {len(parsed_df)}")
    print(f"After keeping latest versions: {len(latest_df)}")

    return latest_df


def create_waterfall_chart(ax, total_items, has_valid_ids_provided, has_valid_ids, mapped_to_kg, title=""):
    """Create a single waterfall chart showing the progression from provided IDs to mapped"""

    # Calculate the differences between each stage
    diff1 = has_valid_ids - has_valid_ids_provided  # Assigned IDs added (usually positive)
    diff2 = has_valid_ids - mapped_to_kg  # Failed to map (usually positive decrease)

    # Define the stages - 2 main columns with 2 labeled difference columns between them
    stages = ['Provided\nIDs', 'Assigned\nIDs', 'Failed to\nmap', 'Mapped\nto KG']
    diffs = [0, diff1, diff2, 0]  # Actual differences
    cumulative = [has_valid_ids_provided, has_valid_ids, has_valid_ids, mapped_to_kg]

    # Create the waterfall chart
    for i, (stage, diff, cum_val) in enumerate(zip(stages, diffs, cumulative)):
        # Convert to percentage of total items for display
        cum_val_pct = (cum_val / total_items) * 100 if total_items > 0 else 0

        if i in [0, 3]:  # Main bars (has provided IDs, mapped to KG)
            if i == 0:
                color = '#808080'  # Grey for has provided IDs
            else:
                color = '#4A90E2'  # Blue for mapped to KG

            ax.bar(i, cum_val_pct, color=color, alpha=0.8, width=0.6)
            ax.text(i, cum_val_pct + 2, f'{cum_val:,}',
                   ha='center', va='bottom', fontweight='bold', fontsize=10)

        elif i in [1, 2]:  # Difference bars (Assigned IDs, Failed to map)
            if i == 1:  # Assigned IDs column
                prev_cum = cumulative[0]  # has_valid_ids_provided
                next_cum = cumulative[1]  # has_valid_ids
            else:  # i == 2, Failed to map column
                prev_cum = cumulative[1]  # has_valid_ids
                next_cum = cumulative[3]  # mapped_to_kg

            prev_cum_pct = (prev_cum / total_items) * 100 if total_items > 0 else 0
            next_cum_pct = (next_cum / total_items) * 100 if total_items > 0 else 0

            # Always show the change, even if it's 0
            diff_pct = abs(diff) / total_items * 100 if total_items > 0 else 0

            if i == 1:  # Assigned IDs - usually an increase
                if diff >= 0:  # Increase or no change
                    color = '#00AA00'  # Green for increases
                    if diff_pct > 0:  # Only draw bar if there's actually a change
                        ax.bar(i, diff_pct, bottom=prev_cum_pct, color=color, alpha=0.8, width=0.6)
                    ax.text(i, max(next_cum_pct, prev_cum_pct) + 2, f'+{diff:,}',
                           ha='center', va='bottom', fontweight='bold', fontsize=9, color='black')
                else:  # Decrease (fewer valid IDs after assignment - shouldn't happen but just in case)
                    color = '#FF0000'  # Red for decreases
                    ax.bar(i, diff_pct, bottom=next_cum_pct, color=color, alpha=0.8, width=0.6)
                    ax.text(i, prev_cum_pct + 2, f'-{abs(diff):,}',
                           ha='center', va='bottom', fontweight='bold', fontsize=9, color='black')
            else:  # i == 2, Failed to map - usually a decrease
                color = '#FF0000'  # Red for decreases (items that failed to map)
                if diff_pct > 0:  # Only draw bar if there's actually a change
                    ax.bar(i, diff_pct, bottom=next_cum_pct, color=color, alpha=0.8, width=0.6)
                ax.text(i, max(prev_cum_pct, next_cum_pct) + 2, f'-{diff:,}',
                       ha='center', va='bottom', fontweight='bold', fontsize=9, color='black')

            # Add connector lines for the difference bars
            ax.plot([i-0.3, i-0.3], [prev_cum_pct, prev_cum_pct], 'k--', alpha=0.5, linewidth=1)
            ax.plot([i-0.3, i+0.3], [prev_cum_pct, prev_cum_pct], 'k--', alpha=0.5, linewidth=1)
            ax.plot([i+0.3, i+0.3], [prev_cum_pct, next_cum_pct], 'k--', alpha=0.5, linewidth=1)

    # Add dotted connecting lines between all columns (like typical waterfall charts)
    heights_pct = [(cumulative[i] / total_items) * 100 if total_items > 0 else 0 for i in range(len(stages))]

    for i in range(len(stages) - 1):
        current_height = heights_pct[i]
        next_height = heights_pct[i+1]

        if i + 1 == 1:  # Connecting to "Assigned IDs" column
            if diffs[1] > 0:  # Increase - line comes in at bottom, goes out at top
                ax.plot([i + 0.3, i + 0.7], [current_height, current_height], 'k:', alpha=0.6, linewidth=1.5, zorder=1)
            else:  # Decrease - line comes in at top, goes out at bottom
                ax.plot([i + 0.3, i + 0.7], [current_height, next_height], 'k:', alpha=0.6, linewidth=1.5, zorder=1)
        elif i + 1 == 2:  # Connecting to "Failed to map" column
            prev_height = heights_pct[1]  # height after assigned IDs
            # This is always a decrease - line comes in at top, goes out at bottom
            ax.plot([i + 0.3, i + 0.7], [prev_height, next_height], 'k:', alpha=0.6, linewidth=1.5, zorder=1)
        elif i + 1 == 3:  # Connecting to "Mapped to KG" column
            prev_height = heights_pct[1]  # height after assigned IDs (before failed to map)
            # Line comes from after failed to map
            ax.plot([i + 0.3, i + 0.7], [next_height, next_height], 'k:', alpha=0.6, linewidth=1.5, zorder=1)

    # Customize the chart
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages, fontsize=9, ha='center')
    ax.set_ylim(0, 110)  # Keep it as percentage scale (0-100% plus some margin)
    ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=10)

    if title:
        ax.set_title(title, fontweight='bold', fontsize=11, pad=10)

    # Add grid for better readability
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)


def create_mapping_type_chart(ax, one_to_one, multi_mappings, title=""):
    """Create a chart showing 1:1 vs multi mappings as percentages"""

    total_mappings = one_to_one + multi_mappings

    if total_mappings == 0:
        ax.text(0.5, 0.5, 'No\nMappings', ha='center', va='center', fontsize=12,
               transform=ax.transAxes, color='gray', fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])
        return

    # Calculate percentages
    one_to_one_pct = (one_to_one / total_mappings) * 100
    multi_pct = (multi_mappings / total_mappings) * 100

    # Create stacked bar
    categories = ['']
    bars1 = ax.bar(categories, one_to_one_pct, color='#A8D8A8', alpha=0.8, label='1:1 Mappings')
    bars2 = ax.bar(categories, multi_pct, bottom=one_to_one_pct, color='#F4C2A8', alpha=0.8, label='Multi Mappings')

    # Add percentage labels
    if one_to_one_pct > 10:  # Only show if segment is large enough
        ax.text(0, one_to_one_pct/2, f'$\\mathbf{{{one_to_one_pct:.1f}\\%}}$\n({one_to_one:,})',
               ha='center', va='center', fontsize=10, fontweight='bold')
    if multi_pct > 10:  # Only show if segment is large enough
        ax.text(0, one_to_one_pct + multi_pct/2, f'$\\mathbf{{{multi_pct:.1f}\\%}}$\n({multi_mappings:,})',
               ha='center', va='center', fontsize=10, fontweight='bold')

    ax.set_ylim(0, 100)
    ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)

    if title:
        ax.set_title(title, fontweight='bold', fontsize=11, pad=10)


def create_id_source_chart(ax, only_provided, only_assigned, both_provided_assigned, title=""):
    """Create a chart showing distribution of ID sources as percentages"""

    total_items = only_provided + only_assigned + both_provided_assigned

    if total_items == 0:
        ax.text(0.5, 0.5, 'No Valid\nIDs', ha='center', va='center', fontsize=12,
               transform=ax.transAxes, color='gray', fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])
        return

    # Calculate percentages
    provided_pct = (only_provided / total_items) * 100
    assigned_pct = (only_assigned / total_items) * 100
    both_pct = (both_provided_assigned / total_items) * 100

    # Create stacked bar
    categories = ['']
    bars1 = ax.bar(categories, provided_pct, color='lightblue', alpha=0.8, label='Only Provided')
    bars2 = ax.bar(categories, both_pct, bottom=provided_pct, color='#D2C3DA', alpha=0.8, label='Both')
    bars3 = ax.bar(categories, assigned_pct, bottom=provided_pct + both_pct, color='#f4cbd5', alpha=0.8, label='Only Assigned')

    # Add percentage labels
    if provided_pct > 8:  # Only show if segment is large enough
        ax.text(0, provided_pct/2, f'$\\mathbf{{{provided_pct:.1f}\\%}}$\n({only_provided:,})',
               ha='center', va='center', fontsize=9)
    if both_pct > 8:  # Only show if segment is large enough
        ax.text(0, provided_pct + both_pct/2, f'$\\mathbf{{{both_pct:.1f}\\%}}$\n({both_provided_assigned:,})',
               ha='center', va='center', fontsize=9)
    if assigned_pct > 8:  # Only show if segment is large enough
        ax.text(0, provided_pct + both_pct + assigned_pct/2, f'$\\mathbf{{{assigned_pct:.1f}\\%}}$\n({only_assigned:,})',
               ha='center', va='center', fontsize=9)

    ax.set_ylim(0, 100)
    ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_axisbelow(True)

    if title:
        ax.set_title(title, fontweight='bold', fontsize=11, pad=10)


def create_waterfall_charts_grid(data_df: pd.DataFrame, kg_name: str, output_dir: Path):
    """Create grid visualization with waterfall charts in each cell"""

    # Get unique datasets and entity types with custom ordering
    dataset_order = ['arivale', 'ukbb', 'israeli10k']
    all_datasets = data_df['dataset'].unique()
    datasets = [d for d in dataset_order if d in all_datasets]

    # Order entity types in reverse alphabetical order
    entity_types = sorted(data_df['entity_type'].unique(), reverse=True)

    # Map names for display
    entity_type_mapping = {
        'clinicallabs': 'Clinical Labs',
        'lipids': 'Lipids',
        'metabolites': 'Metabolites',
        'proteins': 'Proteins'
    }

    datasets_display = []
    for dataset in datasets:
        if dataset.lower() == 'ukbb':
            datasets_display.append('UKBB')
        elif dataset.lower() == 'israeli10k':
            datasets_display.append('HPP')
        else:
            datasets_display.append(dataset.title())

    entity_types_display = [entity_type_mapping.get(et, et.title()) for et in entity_types]

    # Create the grid (entity types on rows, datasets on columns)
    fig, axes = plt.subplots(len(entity_types), len(datasets),
                            figsize=(len(datasets) * 4, len(entity_types) * 3))

    if len(entity_types) == 1:
        axes = axes.reshape(1, -1)
    if len(datasets) == 1:
        axes = axes.reshape(-1, 1)

    # Process each cell
    for i, (entity_type, entity_type_display) in enumerate(zip(entity_types, entity_types_display)):
        for j, (dataset, dataset_display) in enumerate(zip(datasets, datasets_display)):
            ax = axes[i, j]

            # Find data for this combination
            mask = (data_df['dataset'] == dataset) & (data_df['entity_type'] == entity_type)
            if mask.any():
                row = data_df[mask].iloc[0]

                # Extract values
                total_items = row['total_items']
                has_valid_ids_provided = row['has_valid_ids_provided']
                has_valid_ids = row['has_valid_ids']
                mapped_to_kg = row['mapped_to_kg']

                # Create waterfall chart
                create_waterfall_chart(ax, total_items, has_valid_ids_provided, has_valid_ids, mapped_to_kg)

            else:
                # No data for this combination
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=14,
                       transform=ax.transAxes, color='gray', fontweight='bold')
                ax.set_xticks([])
                ax.set_yticks([])

            # Add cell labels
            if i == 0:  # Top row
                ax.set_title(dataset_display, fontsize=12, fontweight='bold', pad=15)
            if j == 0:  # Left column
                ax.set_ylabel(f'{entity_type_display}\n\nPercentage (%)', fontsize=11, fontweight='bold')



    # Adjust layout for better spacing
    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.08)

    # Save the plot
    output_file = output_dir / f'a_waterfall_charts_grid_{kg_name.lower()}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved waterfall charts grid: {output_file}")

    # Also save as PDF
    output_file_pdf = output_dir / f'a_waterfall_charts_grid_{kg_name.lower()}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved waterfall charts grid (PDF): {output_file_pdf}")

    plt.close()


def create_mapping_type_grid(data_df: pd.DataFrame, kg_name: str, output_dir: Path):
    """Create grid visualization showing 1:1 vs multi mapping percentages"""

    # Get unique datasets and entity types with custom ordering
    dataset_order = ['arivale', 'ukbb', 'israeli10k']
    all_datasets = data_df['dataset'].unique()
    datasets = [d for d in dataset_order if d in all_datasets]

    # Order entity types in reverse alphabetical order
    entity_types = sorted(data_df['entity_type'].unique(), reverse=True)

    # Map names for display
    entity_type_mapping = {
        'clinicallabs': 'Clinical Labs',
        'lipids': 'Lipids',
        'metabolites': 'Metabolites',
        'proteins': 'Proteins'
    }

    datasets_display = []
    for dataset in datasets:
        if dataset.lower() == 'ukbb':
            datasets_display.append('UKBB')
        elif dataset.lower() == 'israeli10k':
            datasets_display.append('HPP')
        else:
            datasets_display.append(dataset.title())

    entity_types_display = [entity_type_mapping.get(et, et.title()) for et in entity_types]

    # Create the grid
    fig, axes = plt.subplots(len(entity_types), len(datasets),
                            figsize=(len(datasets) * 3, len(entity_types) * 2.5))

    if len(entity_types) == 1:
        axes = axes.reshape(1, -1)
    if len(datasets) == 1:
        axes = axes.reshape(-1, 1)

    # Process each cell
    for i, (entity_type, entity_type_display) in enumerate(zip(entity_types, entity_types_display)):
        for j, (dataset, dataset_display) in enumerate(zip(datasets, datasets_display)):
            ax = axes[i, j]

            # Find data for this combination
            mask = (data_df['dataset'] == dataset) & (data_df['entity_type'] == entity_type)
            if mask.any():
                row = data_df[mask].iloc[0]

                # Extract mapping values
                one_to_one = row['one_to_one_mappings']
                multi_mappings = row['multi_mappings']

                # Create mapping type chart
                create_mapping_type_chart(ax, one_to_one, multi_mappings)

            else:
                # No data for this combination
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=14,
                       transform=ax.transAxes, color='gray', fontweight='bold')
                ax.set_xticks([])
                ax.set_yticks([])

            # Add cell labels
            if i == 0:  # Top row
                ax.set_title(dataset_display, fontsize=12, fontweight='bold', pad=15)
            if j == 0:  # Left column
                ax.set_ylabel(f'{entity_type_display}\n\nPercentage (%)', fontsize=11, fontweight='bold')


    # Add legend
    legend_elements = [
        patches.Patch(color='#A8D8A8', alpha=0.8, label='1:1 Mappings'),
        patches.Patch(color='#F4C2A8', alpha=0.8, label='Multi Mappings')
    ]
    fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.02),
              ncol=2, fontsize=11)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.08)

    # Save the plot
    output_file = output_dir / f'a_mapping_types_grid_{kg_name.lower()}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved mapping types grid: {output_file}")

    # Also save as PDF
    output_file_pdf = output_dir / f'a_mapping_types_grid_{kg_name.lower()}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved mapping types grid (PDF): {output_file_pdf}")

    plt.close()


def create_id_source_grid(data_df: pd.DataFrame, kg_name: str, output_dir: Path):
    """Create grid visualization showing ID source distribution"""

    # Get unique datasets and entity types with custom ordering
    dataset_order = ['arivale', 'ukbb', 'israeli10k']
    all_datasets = data_df['dataset'].unique()
    datasets = [d for d in dataset_order if d in all_datasets]

    # Order entity types in reverse alphabetical order
    entity_types = sorted(data_df['entity_type'].unique(), reverse=True)

    # Map names for display
    entity_type_mapping = {
        'clinicallabs': 'Clinical Labs',
        'lipids': 'Lipids',
        'metabolites': 'Metabolites',
        'proteins': 'Proteins'
    }

    datasets_display = []
    for dataset in datasets:
        if dataset.lower() == 'ukbb':
            datasets_display.append('UKBB')
        elif dataset.lower() == 'israeli10k':
            datasets_display.append('HPP')
        else:
            datasets_display.append(dataset.title())

    entity_types_display = [entity_type_mapping.get(et, et.title()) for et in entity_types]

    # Create the grid
    fig, axes = plt.subplots(len(entity_types), len(datasets),
                            figsize=(len(datasets) * 3, len(entity_types) * 2.5))

    if len(entity_types) == 1:
        axes = axes.reshape(1, -1)
    if len(datasets) == 1:
        axes = axes.reshape(-1, 1)

    # Process each cell
    for i, (entity_type, entity_type_display) in enumerate(zip(entity_types, entity_types_display)):
        for j, (dataset, dataset_display) in enumerate(zip(datasets, datasets_display)):
            ax = axes[i, j]

            # Find data for this combination
            mask = (data_df['dataset'] == dataset) & (data_df['entity_type'] == entity_type)
            if mask.any():
                row = data_df[mask].iloc[0]

                # Extract ID source values
                only_provided = row['has_only_provided_ids']
                only_assigned = row['has_only_assigned_ids']
                both_provided_assigned = row['has_both_provided_and_assigned_ids']

                # Create ID source chart
                create_id_source_chart(ax, only_provided, only_assigned, both_provided_assigned)

            else:
                # No data for this combination
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=14,
                       transform=ax.transAxes, color='gray', fontweight='bold')
                ax.set_xticks([])
                ax.set_yticks([])

            # Add cell labels
            if i == 0:  # Top row
                ax.set_title(dataset_display, fontsize=12, fontweight='bold', pad=15)
            if j == 0:  # Left column
                ax.set_ylabel(f'{entity_type_display}\n\nPercentage (%)', fontsize=11, fontweight='bold')


    # Add legend
    legend_elements = [
        patches.Patch(color='lightblue', alpha=0.6, label='Only Provided IDs'),
        patches.Patch(color='#f4cbd5', alpha=0.6, label='Only Assigned IDs'),
        patches.Patch(color='white', edgecolor='black', label='Both (Intersection)')
    ]
    fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.02),
              ncol=3, fontsize=10)

    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.90, bottom=0.08)

    # Save the plot
    output_file = output_dir / f'a_id_sources_grid_{kg_name.lower()}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved ID sources grid: {output_file}")

    # Also save as PDF
    output_file_pdf = output_dir / f'a_id_sources_grid_{kg_name.lower()}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved ID sources grid (PDF): {output_file_pdf}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Create waterfall charts grid for mapping pipeline visualization')
    parser.add_argument('kg_name', help='Knowledge graph name (e.g., KRAKEN)')
    parser.add_argument('--input', default='src/mapping/results/kraken/a_summary.tsv',
                       help='Input TSV file path (default: src/mapping/results/kraken/a_summary.tsv)')
    parser.add_argument('--output', default='src/mapping',
                       help='Output directory (default: src/mapping)')
    args = parser.parse_args()

    # File paths
    input_file = Path(args.input)
    output_dir = Path(args.output)

    if not input_file.exists():
        print(f"Error: Input file {input_file} not found!")
        return

    # Load and process data
    print(f"Loading data from: {input_file}")
    data_df = load_and_process_summary_data(input_file)

    print(f"\nData shape: {data_df.shape}")
    print(f"Datasets: {sorted(data_df['dataset'].unique())}")
    print(f"Entity types: {sorted(data_df['entity_type'].unique())}")

    # Create all three visualizations
    print("\nCreating waterfall charts grid...")
    create_waterfall_charts_grid(data_df, args.kg_name, output_dir)

    print("\nCreating mapping types grid...")
    create_mapping_type_grid(data_df, args.kg_name, output_dir)

    print("\nCreating ID sources grid...")
    create_id_source_grid(data_df, args.kg_name, output_dir)

    print(f"\nAll visualizations complete!")
    print(f"\nUsage examples:")
    print(f"  python src/mapping/visualize_waterfall_grid.py KRAKEN")
    print(f"  python src/mapping/visualize_waterfall_grid.py KG2 --input src/mapping/results/kg2/a_summary.tsv")


if __name__ == "__main__":
    main()