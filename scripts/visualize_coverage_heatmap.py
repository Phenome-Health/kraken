#!/usr/bin/env python3
"""
Visualize coverage data as a heatmap with detailed cell annotations.
Shows percentages and mapped/total counts for each entity type and dataset.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse


def load_and_process_data(file_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load TSV data and process into coverage matrix and annotation matrix"""
    
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Remove empty rows
    df = df.dropna(how='all')
    
    # Get unique dataset names (everything before '_coverage', '_mapped', '_total')
    datasets = []
    for col in df.columns:
        if col.endswith('_coverage'):
            dataset = col.replace('_coverage', '')
            datasets.append(dataset)
    
    print(f"Found datasets: {datasets}")
    print(f"Entity types: {df['Entity_Type'].dropna().tolist()}")
    
    # Create coverage matrix and annotation matrix
    entity_types = df['Entity_Type'].dropna().tolist()
    coverage_matrix = pd.DataFrame(index=entity_types, columns=datasets, dtype=float)
    annotation_matrix = pd.DataFrame(index=entity_types, columns=datasets, dtype=str)
    
    for dataset in datasets:
        coverage_col = f"{dataset}_coverage"
        mapped_col = f"{dataset}_mapped"
        total_col = f"{dataset}_total"
        
        for idx, entity_type in enumerate(entity_types):
            if idx < len(df):
                coverage = df.loc[idx, coverage_col]
                mapped = df.loc[idx, mapped_col]
                total = df.loc[idx, total_col]
                
                # Handle N/A values
                if pd.isna(coverage) or coverage == 'N/A':
                    coverage_matrix.loc[entity_type, dataset] = np.nan
                    annotation_matrix.loc[entity_type, dataset] = 'N/A'
                else:
                    try:
                        coverage_val = float(coverage)
                        mapped_val = int(mapped) if not pd.isna(mapped) else 0
                        total_val = int(total) if not pd.isna(total) else 0
                        
                        coverage_matrix.loc[entity_type, dataset] = coverage_val
                        annotation_matrix.loc[entity_type, dataset] = f"{coverage_val:.1f}%\n{mapped_val:,} / {total_val:,}"
                    except (ValueError, TypeError):
                        coverage_matrix.loc[entity_type, dataset] = np.nan
                        annotation_matrix.loc[entity_type, dataset] = 'N/A'
    
    return coverage_matrix, annotation_matrix


def create_coverage_heatmap(coverage_matrix: pd.DataFrame, annotation_matrix: pd.DataFrame, output_dir: Path):
    """Create and save the coverage heatmap"""
    
    # Set up the plot
    plt.figure(figsize=(12, 8))
    
    # Create heatmap with custom colormap
    mask = coverage_matrix.isna()
    
    # Use a colormap that goes from red (low coverage) to green (high coverage)
    cmap = sns.diverging_palette(10, 130, as_cmap=True)  # Red to green
    
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
    
    # Customize the plot
    plt.title('Entity Type Coverage Across Datasets', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Dataset', fontsize=12, fontweight='bold')
    plt.ylabel('Entity Type', fontsize=12, fontweight='bold')
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    
    # Save the plot
    output_file = output_dir / 'coverage_heatmap.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved coverage heatmap: {output_file}")
    
    # Also save as PDF for vector graphics
    output_file_pdf = output_dir / 'coverage_heatmap.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Saved coverage heatmap (PDF): {output_file_pdf}")
    
    plt.close()


def generate_summary_stats(coverage_matrix: pd.DataFrame, output_dir: Path):
    """Generate and save summary statistics"""
    
    print("\n=== COVERAGE SUMMARY STATISTICS ===")
    
    # Overall statistics
    valid_coverage = coverage_matrix.dropna()
    if not valid_coverage.empty:
        mean_coverage = valid_coverage.values.flatten().mean()
        print(f"Average coverage across all datasets and entity types: {mean_coverage:.1f}%")
        
        # Dataset-wise statistics
        print(f"\nDataset-wise average coverage:")
        for dataset in coverage_matrix.columns:
            dataset_coverage = coverage_matrix[dataset].dropna()
            if not dataset_coverage.empty:
                avg = dataset_coverage.mean()
                print(f"  {dataset}: {avg:.1f}%")
        
        # Entity type-wise statistics
        print(f"\nEntity type-wise average coverage:")
        for entity_type in coverage_matrix.index:
            entity_coverage = coverage_matrix.loc[entity_type].dropna()
            if not entity_coverage.empty:
                avg = entity_coverage.mean()
                print(f"  {entity_type}: {avg:.1f}%")
    
    # Save summary to file
    summary_file = output_dir / 'coverage_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("=== COVERAGE SUMMARY STATISTICS ===\n\n")
        
        if not valid_coverage.empty:
            f.write(f"Average coverage across all datasets and entity types: {mean_coverage:.1f}%\n\n")
            
            f.write("Dataset-wise average coverage:\n")
            for dataset in coverage_matrix.columns:
                dataset_coverage = coverage_matrix[dataset].dropna()
                if not dataset_coverage.empty:
                    avg = dataset_coverage.mean()
                    f.write(f"  {dataset}: {avg:.1f}%\n")
            
            f.write(f"\nEntity type-wise average coverage:\n")
            for entity_type in coverage_matrix.index:
                entity_coverage = coverage_matrix.loc[entity_type].dropna()
                if not entity_coverage.empty:
                    avg = entity_coverage.mean()
                    f.write(f"  {entity_type}: {avg:.1f}%\n")
    
    print(f"Saved summary: {summary_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize coverage data as a heatmap')
    parser.add_argument('--input', default='scripts/coverage_data.tsv', 
                       help='Input TSV file path (default: scripts/coverage_data.tsv)')
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
    create_coverage_heatmap(coverage_matrix, annotation_matrix, output_dir)
    
    # Generate summary statistics
    generate_summary_stats(coverage_matrix, output_dir)
    
    print(f"\nUsage tips:")
    print(f"  python scripts/visualize_coverage_heatmap.py                    # Use default files")
    print(f"  python scripts/visualize_coverage_heatmap.py --input data.tsv   # Use custom input file")


if __name__ == "__main__":
    main()