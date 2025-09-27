#!/usr/bin/env python3
"""
Visualize mapping summary data as a heatmap with detailed cell annotations.
Shows percentages and mapped/total counts for each entity type and dataset.
Processes data from mapping results summary files.
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


def load_and_process_data(file_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load TSV data and process into coverage matrix and annotation matrix"""

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

    print(f"Original entries: {len(parsed_df)}")
    print(f"After keeping latest versions: {len(latest_df)}")

    # Get unique datasets and entity types with custom ordering
    # Order datasets: arivale, ukbb, israeli10k
    all_datasets = latest_df['dataset'].unique()
    dataset_order = ['arivale', 'ukbb', 'israeli10k']
    datasets = [d for d in dataset_order if d in all_datasets]

    # Order entity types in reverse alphabetical order
    entity_types = sorted(latest_df['entity_type'].unique(), reverse=True)

    # Map entity type names for display
    entity_type_mapping = {
        'clinicallabs': 'Clinical Labs',
        'lipids': 'Lipids',
        'metabolites': 'Metabolites',
        'proteins': 'Proteins'
    }

    entity_types_display = [entity_type_mapping.get(et, et.title()) for et in entity_types]

    print(f"Datasets: {datasets}")
    print(f"Entity types: {entity_types}")
    print(f"Entity types display: {entity_types_display}")

    # Create coverage matrix and annotation matrix
    coverage_matrix = pd.DataFrame(index=entity_types_display, columns=datasets, dtype=float)
    annotation_matrix = pd.DataFrame(index=entity_types_display, columns=datasets, dtype=str)

    # Fill matrices
    for _, row in latest_df.iterrows():
        dataset = row['dataset']
        entity_type = row['entity_type']
        entity_type_display = entity_type_mapping.get(entity_type, entity_type.title())

        total_items = row['total_items']
        mapped_to_kg = row['mapped_to_kg']

        # Calculate percentage
        if total_items > 0:
            coverage_percentage = (mapped_to_kg / total_items) * 100
            coverage_matrix.loc[entity_type_display, dataset] = coverage_percentage
            annotation_matrix.loc[entity_type_display, dataset] = f"{coverage_percentage:.1f}%\n\n{mapped_to_kg:,} / {total_items:,}"
        else:
            coverage_matrix.loc[entity_type_display, dataset] = -1  # Use -1 as placeholder for N/A
            annotation_matrix.loc[entity_type_display, dataset] = 'N/A'

    # Fill missing combinations with N/A
    for entity_type_display in entity_types_display:
        for dataset in datasets:
            if pd.isna(coverage_matrix.loc[entity_type_display, dataset]):
                coverage_matrix.loc[entity_type_display, dataset] = -1
                annotation_matrix.loc[entity_type_display, dataset] = 'N/A'

    return coverage_matrix, annotation_matrix


def create_coverage_heatmap(coverage_matrix: pd.DataFrame, annotation_matrix: pd.DataFrame, output_dir: Path, kg_name: str):
    """Create and save the coverage heatmap"""

    # Set up the plot
    plt.figure(figsize=(12, 8))

    # Create heatmap with custom colormap
    # Use a colormap that goes from red (low coverage) to green (high coverage)
    cmap = sns.diverging_palette(10, 130, as_cmap=True)  # Red to green

    # Create mask for N/A values (where coverage_matrix == -1)
    mask = (coverage_matrix == -1)

    ax = sns.heatmap(
        coverage_matrix,
        annot=annotation_matrix,
        fmt='',
        cmap=cmap,
        center=50,  # Center the colormap at 50%
        vmin=0,
        vmax=100,
        mask=mask,
        cbar_kws={'label': 'Coverage Percentage'},
        annot_kws={'size': 10, 'ha': 'center', 'va': 'center'},
        square=True,
        linewidths=0.5,
        linecolor='white'
    )

    # Manually add N/A annotations for masked cells
    for i, row_name in enumerate(coverage_matrix.index):
        for j, col_name in enumerate(coverage_matrix.columns):
            if coverage_matrix.loc[row_name, col_name] == -1:
                ax.text(j + 0.5, i + 0.5, 'N/A',
                       ha='center', va='center', fontsize=10,
                       color='black')

    # Customize the plot
    plt.title(f'Dataset-to-{kg_name} Mapping Summary - Coverage', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Dataset', fontsize=12, fontweight='bold')
    plt.ylabel('Entity Type', fontsize=12, fontweight='bold')

    # Keep x-axis labels horizontal
    plt.xticks(rotation=0, ha='center')
    plt.yticks(rotation=0)

    # Add footnote explaining the data
    plt.figtext(0.6, -0.03, f'Values show percentage of items successfully mapped to {kg_name}.',
                fontsize=9, style='italic', ha='center', va='bottom')

    # Adjust layout to prevent label cutoff and accommodate footnote
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.08)  # Make room for footnote

    # Save the plot
    output_file = output_dir / f'a_summary_coverage_heatmap_{kg_name.lower()}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved summary coverage heatmap: {output_file}")

    # Also save as PDF for vector graphics
    output_file_pdf = output_dir / f'a_summary_coverage_heatmap_{kg_name.lower()}.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved summary coverage heatmap (PDF): {output_file_pdf}")

    plt.close()


def generate_summary_stats(coverage_matrix: pd.DataFrame, output_dir: Path, kg_name: str):
    """Generate and save summary statistics"""

    print("\n=== SUMMARY COVERAGE STATISTICS ===")

    # Overall statistics (excluding N/A values)
    valid_coverage = coverage_matrix[coverage_matrix != -1]
    if not valid_coverage.empty:
        mean_coverage = valid_coverage.values.flatten().mean()
        print(f"Average coverage across all datasets and entity types: {mean_coverage:.1f}%")

        # Dataset-wise statistics
        print(f"\nDataset-wise average coverage:")
        for dataset in coverage_matrix.columns:
            dataset_coverage = coverage_matrix[dataset][coverage_matrix[dataset] != -1]
            if not dataset_coverage.empty:
                avg = dataset_coverage.mean()
                print(f"  {dataset}: {avg:.1f}%")

        # Entity type-wise statistics
        print(f"\nEntity type-wise average coverage:")
        for entity_type in coverage_matrix.index:
            entity_coverage = coverage_matrix.loc[entity_type][coverage_matrix.loc[entity_type] != -1]
            if not entity_coverage.empty:
                avg = entity_coverage.mean()
                print(f"  {entity_type}: {avg:.1f}%")

    # Save summary to file
    summary_file = output_dir / f'a_summary_coverage_stats_{kg_name.lower()}.txt'
    with open(summary_file, 'w') as f:
        f.write("=== SUMMARY COVERAGE STATISTICS ===\n\n")

        if not valid_coverage.empty:
            f.write(f"Average coverage across all datasets and entity types: {mean_coverage:.1f}%\n\n")

            f.write("Dataset-wise average coverage:\n")
            for dataset in coverage_matrix.columns:
                dataset_coverage = coverage_matrix[dataset][coverage_matrix[dataset] != -1]
                if not dataset_coverage.empty:
                    avg = dataset_coverage.mean()
                    f.write(f"  {dataset}: {avg:.1f}%\n")

            f.write(f"\nEntity type-wise average coverage:\n")
            for entity_type in coverage_matrix.index:
                entity_coverage = coverage_matrix.loc[entity_type][coverage_matrix.loc[entity_type] != -1]
                if not entity_coverage.empty:
                    avg = entity_coverage.mean()
                    f.write(f"  {entity_type}: {avg:.1f}%\n")

    print(f"Saved summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize mapping summary data as a heatmap')
    parser.add_argument('kg_name', help='Knowledge graph name (e.g., KRAKEN)')
    parser.add_argument('--input', default='src/mapping/results/kraken/a_summary.tsv',
                       help='Input TSV file path (default: src/mapping/results/kraken/a_summary.tsv)')
    parser.add_argument('--output', default='scripts',
                       help='Output directory (default: scripts)')
    args = parser.parse_args()

    # File paths
    input_file = Path(args.input)
    output_dir = Path(args.output)

    if not input_file.exists():
        print(f"Error: Input file {input_file} not found!")
        return

    # Load and process data
    print(f"Loading data from: {input_file}")
    coverage_matrix, annotation_matrix = load_and_process_data(input_file)

    print(f"\nCoverage matrix shape: {coverage_matrix.shape}")
    print(f"Datasets: {list(coverage_matrix.columns)}")
    print(f"Entity types: {list(coverage_matrix.index)}")

    # Create visualizations
    create_coverage_heatmap(coverage_matrix, annotation_matrix, output_dir, args.kg_name)

    # Generate summary statistics
    generate_summary_stats(coverage_matrix, output_dir, args.kg_name)


if __name__ == "__main__":
    main()