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
                'has_valid_ids_provided': row['has_valid_ids_provided'],
                'mapped_to_kg': row['mapped_to_kg'],
                'one_to_one_mappings': row['one_to_one_mappings'],
                'one_to_many_mappings': row['one_to_many_mappings']
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
        else:
            datasets_display.append(dataset.title())

    entity_types_display = [entity_type_mapping.get(et, et.title()) for et in entity_types]

    # Create the grid
    fig, axes = plt.subplots(len(datasets), len(entity_types), figsize=(len(entity_types) * 3, len(datasets) * 2))
    if len(datasets) == 1:
        axes = axes.reshape(1, -1)
    if len(entity_types) == 1:
        axes = axes.reshape(-1, 1)

    # Process each cell
    for i, (dataset, dataset_display) in enumerate(zip(datasets, datasets_display)):
        for j, (entity_type, entity_type_display) in enumerate(zip(entity_types, entity_types_display)):
            ax = axes[i, j]

            # Find data for this combination
            mask = (data_df['dataset'] == dataset) & (data_df['entity_type'] == entity_type)
            if mask.any():
                row = data_df[mask].iloc[0]

                # Extract values
                total_items = row['total_items']
                has_valid_ids = row['has_valid_ids']
                has_valid_ids_provided = row['has_valid_ids_provided']
                mapped_to_kg = row['mapped_to_kg']
                one_to_one = row['one_to_one_mappings']
                one_to_many = row['one_to_many_mappings']

                # Calculate percentages
                valid_ids_pct = (has_valid_ids / total_items) * 100 if total_items > 0 else 0
                provided_ids_pct = (has_valid_ids_provided / has_valid_ids) * 100 if has_valid_ids > 0 else 0
                mapped_pct = (mapped_to_kg / total_items) * 100 if total_items > 0 else 0
                one_to_one_pct = (one_to_one / mapped_to_kg) * 100 if mapped_to_kg > 0 else 0
                one_to_many_pct = (one_to_many / mapped_to_kg) * 100 if mapped_to_kg > 0 else 0

                # Create three bars
                bar_width = 0.25
                x_positions = [0.1, 0.4, 0.7]

                # Bar 1: Total items baseline (always 100%)
                ax.bar(x_positions[0], 100, bar_width, color='lightgray', alpha=0.7, label='Total Items')

                # Bar 2: Valid IDs breakdown
                # Bottom part: has_valid_ids_provided (light blue)
                if has_valid_ids > 0:
                    provided_portion = (provided_ids_pct / 100) * valid_ids_pct
                    ax.bar(x_positions[1], provided_portion, bar_width, color='lightblue', alpha=0.8)
                    # Top part: remaining valid IDs (pale purple)
                    remaining_portion = valid_ids_pct - provided_portion
                    ax.bar(x_positions[1], remaining_portion, bar_width, color='thistle', alpha=0.8,
                          bottom=provided_portion)

                    # Add count labels inside the stacked parts
                    if provided_portion > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[1], provided_portion/2, f'{has_valid_ids_provided:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                    if remaining_portion > 8:  # Only show if chunk is big enough
                        assigned_ids = has_valid_ids - has_valid_ids_provided
                        ax.text(x_positions[1], provided_portion + remaining_portion/2, f'{assigned_ids:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                else:
                    ax.bar(x_positions[1], valid_ids_pct, bar_width, color='thistle', alpha=0.8)

                # Bar 3: Mapping breakdown
                # Bottom part: one_to_one (light green)
                if mapped_to_kg > 0:
                    one_to_one_portion = (one_to_one_pct / 100) * mapped_pct
                    one_to_many_portion = (one_to_many_pct / 100) * mapped_pct

                    ax.bar(x_positions[2], one_to_one_portion, bar_width, color='lightgreen', alpha=0.8)
                    # Top part: one_to_many (yellow)
                    ax.bar(x_positions[2], one_to_many_portion, bar_width, color='khaki', alpha=0.8,
                          bottom=one_to_one_portion)

                    # Add count labels inside the stacked parts
                    if one_to_one_portion > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[2], one_to_one_portion/2, f'{one_to_one:,}',
                               ha='center', va='center', fontsize=7, fontweight='bold', color='black')
                    if one_to_many_portion > 8:  # Only show if chunk is big enough
                        ax.text(x_positions[2], one_to_one_portion + one_to_many_portion/2, f'{one_to_many:,}',
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

            # Add cell label (dataset - entity type)
            if i == 0:  # Top row
                ax.set_title(entity_type_display, fontsize=10, fontweight='bold', pad=10)
            if j == 0:  # Left column
                ax.set_ylabel(dataset_display, fontsize=10, fontweight='bold')

    # Add overall title with more space
    fig.suptitle(f'{kg_name} Mapping Pipeline Breakdown', fontsize=16, fontweight='bold', y=0.98)

    # Add legend
    legend_elements = [
        patches.Patch(color='lightgray', alpha=0.7, label='Total Items'),
        patches.Patch(color='lightblue', alpha=0.8, label='Provided IDs'),
        patches.Patch(color='thistle', alpha=0.8, label='Assigned IDs'),
        patches.Patch(color='lightgreen', alpha=0.8, label='1:1 KG Mappings'),
        patches.Patch(color='khaki', alpha=0.8, label='1:Many KG Mappings')
    ]
    fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, 0.02), ncol=5, fontsize=9)

    # Adjust layout to give title more space
    plt.tight_layout()
    plt.subplots_adjust(top=0.85, bottom=0.15)

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