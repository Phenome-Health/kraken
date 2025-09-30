#!/usr/bin/env python3
"""
Visualize mapping summary data as a grid of stacked bar charts.
Each cell contains mini stacked bar charts showing ID validation and mapping breakdown.
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
                'has_only_provided_ids': row['has_only_provided_ids'],
                'has_both_provided_and_assigned_ids': row['has_both_provided_and_assigned_ids'],
                'has_only_assigned_ids': row['has_only_assigned_ids'],
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


def create_stacked_bars_grid(data_df: pd.DataFrame, kg_name: str, output_dir: Path):
    """Create grid visualization with stacked bar charts in each cell"""

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

    # Create the grid (swapped: entity types on rows, datasets on columns)
    fig, axes = plt.subplots(len(entity_types), len(datasets), figsize=(len(datasets) * 3, len(entity_types) * 2))
    if len(entity_types) == 1:
        axes = axes.reshape(1, -1)
    if len(datasets) == 1:
        axes = axes.reshape(-1, 1)

    # Process each cell (swapped indices: entity types on rows, datasets on columns)
    for i, (entity_type, entity_type_display) in enumerate(zip(entity_types, entity_types_display)):
        for j, (dataset, dataset_display) in enumerate(zip(datasets, datasets_display)):
            ax = axes[i, j]

            # Find data for this combination
            mask = (data_df['dataset'] == dataset) & (data_df['entity_type'] == entity_type)
            if mask.any():
                row = data_df[mask].iloc[0]

                # Extract values
                total_items = row['total_items']
                has_valid_ids = row['has_valid_ids']
                has_only_provided_ids = row['has_only_provided_ids']
                has_both_provided_and_assigned_ids = row['has_both_provided_and_assigned_ids']
                has_only_assigned_ids = row['has_only_assigned_ids']
                mapped_to_kg = row['mapped_to_kg']
                one_to_one = row['one_to_one_mappings']
                multi = row['multi_mappings']

                # Calculate percentages
                valid_ids_pct = (has_valid_ids / total_items) * 100 if total_items > 0 else 0
                only_provided_pct = (has_only_provided_ids / total_items) * 100 if total_items > 0 else 0
                both_pct = (has_both_provided_and_assigned_ids / total_items) * 100 if total_items > 0 else 0
                only_assigned_pct = (has_only_assigned_ids / total_items) * 100 if total_items > 0 else 0
                mapped_pct = (mapped_to_kg / total_items) * 100 if total_items > 0 else 0
                one_to_one_pct = (one_to_one / mapped_to_kg) * 100 if mapped_to_kg > 0 else 0
                multi_pct = (multi / mapped_to_kg) * 100 if mapped_to_kg > 0 else 0

                # Create three bars
                bar_width = 0.25
                x_positions = [0.1, 0.4, 0.7]

                # Bar 1: Total items baseline (always 100%)
                ax.bar(x_positions[0], 100, bar_width, color='lightgray', alpha=0.7, label='Total Items')

                # Bar 2: Valid IDs breakdown with three chunks
                if has_valid_ids > 0:
                    # Bottom chunk: only provided IDs (light blue)
                    ax.bar(x_positions[1], only_provided_pct, bar_width, color='lightblue', alpha=0.8)
                    # Middle chunk: both provided and assigned IDs (muted purple)
                    ax.bar(x_positions[1], both_pct, bar_width, color='#D2C3DA', alpha=0.8,
                          bottom=only_provided_pct)
                    # Top chunk: only assigned IDs (muted pink)
                    ax.bar(x_positions[1], only_assigned_pct, bar_width, color='#f4cbd5', alpha=0.8,
                          bottom=only_provided_pct + both_pct)

                    # Add count labels inside the stacked parts
                    if only_provided_pct > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[1], only_provided_pct/2, f'{has_only_provided_ids:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                    if both_pct > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[1], only_provided_pct + both_pct/2, f'{has_both_provided_and_assigned_ids:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                    if only_assigned_pct > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[1], only_provided_pct + both_pct + only_assigned_pct/2, f'{has_only_assigned_ids:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                else:
                    ax.bar(x_positions[1], valid_ids_pct, bar_width, color='lightblue', alpha=0.8)

                # Bar 3: Mapping breakdown with two chunks
                if mapped_to_kg > 0:
                    one_to_one_portion = (one_to_one_pct / 100) * mapped_pct
                    multi_portion = (multi_pct / 100) * mapped_pct

                    # Bottom chunk: one_to_one (light green)
                    ax.bar(x_positions[2], one_to_one_portion, bar_width, color='#a6f1a6', alpha=0.8)
                    # Top chunk: multi mappings (khaki)
                    ax.bar(x_positions[2], multi_portion, bar_width, color='khaki', alpha=0.8,
                          bottom=one_to_one_portion)

                    # Add count labels inside the stacked parts
                    if one_to_one_portion > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[2], one_to_one_portion/2, f'{one_to_one:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                    if multi_portion > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[2], one_to_one_portion + multi_portion/2, f'{multi:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')

                # Add count labels below each bar
                ax.text(x_positions[0], -10, f'{total_items:,}', ha='center', va='top', fontsize=8, rotation=0)
                ax.text(x_positions[1], -10, f'{has_valid_ids:,}', ha='center', va='top', fontsize=8, rotation=0)
                ax.text(x_positions[2], -10, f'{mapped_to_kg:,}', ha='center', va='top', fontsize=8, rotation=0)

            else:
                # No data for this combination
                ax.text(0.5, 0.5, 'N/A', ha='center', va='center', fontsize=12,
                       transform=ax.transAxes, color='gray')

            # Customize each subplot
            ax.set_ylim(0, 110)
            ax.set_xlim(0, 1)
            ax.set_xticks([])
            ax.set_yticks([0, 25, 50, 75, 100])
            ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=8)

            # Add grid
            ax.grid(True, alpha=0.3, axis='y')

            # Add cell label (entity type - dataset)
            if i == 0:  # Top row
                ax.set_title(dataset_display, fontsize=10, fontweight='bold', pad=10)
            if j == 0:  # Left column
                ax.set_ylabel(entity_type_display, fontsize=10, fontweight='bold')

    # Add overall title with less space
    fig.suptitle(f'{kg_name} Mapping Results Breakdown', fontsize=16, fontweight='bold', y=0.96)

    # Add legend
    legend_elements = [
        patches.Patch(color='lightgray', alpha=0.7, label='Total Items'),
        patches.Patch(color='lightblue', alpha=0.8, label='Only Provided IDs'),
        patches.Patch(color='#D2C3DA', alpha=0.8, label='Both Provided & Assigned IDs'),
        patches.Patch(color='#f4cbd5', alpha=0.8, label='Only Assigned IDs'),
        patches.Patch(color='#a6f1a6', alpha=0.8, label='1:1 KG Mappings'),
        patches.Patch(color='khaki', alpha=0.8, label='Multi KG Mappings')
    ]
    fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.05), ncol=3, fontsize=9)

    # Adjust layout for better spacing
    plt.tight_layout()
    plt.subplots_adjust(top=0.88, bottom=0.12)

    # Save the plot
    output_file = output_dir / f'a_stacked_bars_grid_{kg_name.lower()}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved stacked bars grid: {output_file}")

    # Also save as PDF
    output_file_pdf = output_dir / f'a_stacked_bars_grid_{kg_name.lower()}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved stacked bars grid (PDF): {output_file_pdf}")

    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Create stacked bar charts grid for mapping pipeline visualization')
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

    # Create visualization
    create_stacked_bars_grid(data_df, args.kg_name, output_dir)

    print(f"\nUsage examples:")
    print(f"  python src/mapping/visualize_stacked_bars_grid.py KRAKEN")
    print(f"  python src/mapping/visualize_stacked_bars_grid.py KG2 --input src/mapping/results/kg2/a_summary.tsv")


if __name__ == "__main__":
    main()