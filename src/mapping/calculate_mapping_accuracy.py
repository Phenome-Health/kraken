#!/usr/bin/env python3
"""
Calculate mapping accuracy by confidence category for full results data.
Analyzes how well assigned IDs match provided IDs across different confidence levels.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import ast


def safe_eval_list(value):
    """Safely evaluate string representation of list"""
    if pd.isna(value) or value == '' or value == '[]':
        return []
    try:
        # Handle the case where it might already be a list
        if isinstance(value, list):
            return value
        # Try to evaluate as Python literal
        return ast.literal_eval(value)
    except (ValueError, SyntaxError):
        # If it fails, try splitting by comma (fallback)
        return [item.strip() for item in str(value).split(',') if item.strip()]


def calculate_set_overlap_accuracy(provided_ids, assigned_ids):
    """Calculate if there's any overlap between provided and assigned ID sets"""
    if not provided_ids or not assigned_ids:
        return None  # Can't calculate accuracy if either is empty

    provided_set = set(provided_ids)
    assigned_set = set(assigned_ids)

    # Return True if there's any overlap
    return len(provided_set.intersection(assigned_set)) > 0


def calculate_majority_id_accuracy(provided_majority, assigned_majority):
    """Calculate if majority canonical IDs match exactly"""
    if pd.isna(provided_majority) or pd.isna(assigned_majority) or \
       provided_majority == '' or assigned_majority == '':
        return None  # Can't calculate accuracy if either is empty

    return provided_majority == assigned_majority


def analyze_accuracy_by_confidence(df):
    """Analyze accuracy metrics by confidence category"""

    # Parse list columns
    print("Parsing ID lists...")
    df['kg_canonical_ids_provided_parsed'] = df['kg_canonical_ids_provided'].apply(safe_eval_list)
    df['kg_canonical_ids_assigned_parsed'] = df['kg_canonical_ids_assigned'].apply(safe_eval_list)

    # Calculate accuracy metrics
    print("Calculating accuracy metrics...")
    df['set_overlap_accuracy'] = df.apply(
        lambda row: calculate_set_overlap_accuracy(
            row['kg_canonical_ids_provided_parsed'],
            row['kg_canonical_ids_assigned_parsed']
        ), axis=1
    )

    df['majority_id_accuracy'] = df.apply(
        lambda row: calculate_majority_id_accuracy(
            row['kg_majority_canonical_id_provided'],
            row['kg_majority_canonical_id_assigned']
        ), axis=1
    )

    # Calculate overall accuracy (across all confidence levels)
    overall_stats = calculate_overall_accuracy(df)

    # Group by confidence and calculate statistics
    confidence_stats = []

    for confidence in ['HIGH', 'MEDIUM', 'LOW']:
        conf_df = df[df['confidence'] == confidence]

        # Set overlap accuracy
        set_overlap_valid = conf_df['set_overlap_accuracy'].notna()
        set_overlap_accurate = conf_df['set_overlap_accuracy'] == True

        set_overlap_total = set_overlap_valid.sum()
        set_overlap_correct = set_overlap_accurate.sum()
        set_overlap_accuracy_rate = (set_overlap_correct / set_overlap_total * 100) if set_overlap_total > 0 else 0

        # Majority ID accuracy
        majority_id_valid = conf_df['majority_id_accuracy'].notna()
        majority_id_accurate = conf_df['majority_id_accuracy'] == True

        majority_id_total = majority_id_valid.sum()
        majority_id_correct = majority_id_accurate.sum()
        majority_id_accuracy_rate = (majority_id_correct / majority_id_total * 100) if majority_id_total > 0 else 0

        confidence_stats.append({
            'confidence': confidence,
            'total_items': len(conf_df),
            'set_overlap_evaluable': set_overlap_total,
            'set_overlap_accurate': set_overlap_correct,
            'set_overlap_accuracy_rate': set_overlap_accuracy_rate,
            'majority_id_evaluable': majority_id_total,
            'majority_id_accurate': majority_id_correct,
            'majority_id_accuracy_rate': majority_id_accuracy_rate
        })

    return pd.DataFrame(confidence_stats), overall_stats


def calculate_overall_accuracy(df):
    """Calculate overall accuracy across all confidence levels"""
    # Set overlap accuracy
    set_overlap_valid = df['set_overlap_accuracy'].notna()
    set_overlap_accurate = df['set_overlap_accuracy'] == True

    set_overlap_total = set_overlap_valid.sum()
    set_overlap_correct = set_overlap_accurate.sum()
    set_overlap_accuracy_rate = (set_overlap_correct / set_overlap_total * 100) if set_overlap_total > 0 else 0

    # Majority ID accuracy
    majority_id_valid = df['majority_id_accuracy'].notna()
    majority_id_accurate = df['majority_id_accuracy'] == True

    majority_id_total = majority_id_valid.sum()
    majority_id_correct = majority_id_accurate.sum()
    majority_id_accuracy_rate = (majority_id_correct / majority_id_total * 100) if majority_id_total > 0 else 0

    return {
        'total_items': len(df),
        'set_overlap_evaluable': set_overlap_total,
        'set_overlap_accurate': set_overlap_correct,
        'set_overlap_accuracy_rate': set_overlap_accuracy_rate,
        'majority_id_evaluable': majority_id_total,
        'majority_id_accurate': majority_id_correct,
        'majority_id_accuracy_rate': majority_id_accuracy_rate
    }


def generate_detailed_breakdown(df):
    """Generate detailed breakdown of accuracy patterns"""
    breakdown = []

    for confidence in ['HIGH', 'MEDIUM', 'LOW']:
        conf_df = df[df['confidence'] == confidence]

        # Count different accuracy patterns
        both_accurate = ((conf_df['set_overlap_accuracy'] == True) &
                        (conf_df['majority_id_accuracy'] == True)).sum()

        set_accurate_only = ((conf_df['set_overlap_accuracy'] == True) &
                           (conf_df['majority_id_accuracy'] == False)).sum()

        majority_accurate_only = ((conf_df['set_overlap_accuracy'] == False) &
                                (conf_df['majority_id_accuracy'] == True)).sum()

        neither_accurate = ((conf_df['set_overlap_accuracy'] == False) &
                          (conf_df['majority_id_accuracy'] == False)).sum()

        # Count cases where only one metric is evaluable
        set_only_evaluable = ((conf_df['set_overlap_accuracy'].notna()) &
                             (conf_df['majority_id_accuracy'].isna())).sum()

        majority_only_evaluable = ((conf_df['set_overlap_accuracy'].isna()) &
                                  (conf_df['majority_id_accuracy'].notna())).sum()

        breakdown.append({
            'confidence': confidence,
            'both_accurate': both_accurate,
            'set_overlap_accurate_only': set_accurate_only,
            'majority_id_accurate_only': majority_accurate_only,
            'neither_accurate': neither_accurate,
            'set_overlap_only_evaluable': set_only_evaluable,
            'majority_id_only_evaluable': majority_only_evaluable
        })

    return pd.DataFrame(breakdown)


def create_accuracy_visualization(stats_df, overall_stats, output_viz_path):
    """Create visualization of accuracy results"""

    # Set up the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot 1: Accuracy rates by confidence category
    confidence_order = ['HIGH', 'MEDIUM', 'LOW']
    stats_ordered = stats_df.set_index('confidence').reindex(confidence_order)

    x_pos = range(len(confidence_order))
    width = 0.35

    # Create bars
    bars1 = ax1.bar([x - width/2 for x in x_pos],
                   stats_ordered['set_overlap_accuracy_rate'],
                   width, label='Set Overlap Accuracy',
                   color='lightblue', alpha=0.8)

    bars2 = ax1.bar([x + width/2 for x in x_pos],
                   stats_ordered['majority_id_accuracy_rate'],
                   width, label='Majority ID Accuracy',
                   color='#FFDAB9', alpha=0.8)

    # Add value labels on bars with counts
    for i, bar in enumerate(bars1):
        height = bar.get_height()
        confidence = confidence_order[i]
        accurate = stats_ordered.loc[confidence, 'set_overlap_accurate']
        evaluable = stats_ordered.loc[confidence, 'set_overlap_evaluable']
        ax1.annotate(f'$\\mathbf{{{height:.1f}\\%}}$\n({accurate:,} / {evaluable:,})',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=11)

    for i, bar in enumerate(bars2):
        height = bar.get_height()
        confidence = confidence_order[i]
        accurate = stats_ordered.loc[confidence, 'majority_id_accurate']
        evaluable = stats_ordered.loc[confidence, 'majority_id_evaluable']
        ax1.annotate(f'$\\mathbf{{{height:.1f}\\%}}$\n({accurate:,} / {evaluable:,})',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=11)

    ax1.set_xlabel('Confidence Category', fontweight='bold')
    ax1.set_ylabel('Accuracy Rate (%)', fontweight='bold')
    ax1.set_title('Mapping Accuracy by Confidence Category', fontweight='bold', fontsize=14)
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels(confidence_order)
    ax1.legend()
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3, axis='y')

    # Plot 2: Overall accuracy comparison
    metrics = ['Set Overlap', 'Majority ID']
    overall_rates = [overall_stats['set_overlap_accuracy_rate'],
                    overall_stats['majority_id_accuracy_rate']]
    colors = ['lightblue', '#FFDAB9']

    bars = ax2.bar(metrics, overall_rates, color=colors, alpha=0.8)

    # Add value labels with counts
    for i, bar in enumerate(bars):
        height = bar.get_height()
        if i == 0:  # Set Overlap
            accurate = overall_stats['set_overlap_accurate']
            evaluable = overall_stats['set_overlap_evaluable']
        else:  # Majority ID
            accurate = overall_stats['majority_id_accurate']
            evaluable = overall_stats['majority_id_evaluable']

        ax2.annotate(f'$\\mathbf{{{height:.1f}\\%}}$\n({accurate:,} / {evaluable:,})',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=12)

    ax2.set_ylabel('Overall Accuracy Rate (%)', fontweight='bold')
    ax2.set_title('Overall Mapping Accuracy', fontweight='bold', fontsize=14)
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3, axis='y')

    # Add sample sizes as text
    ax2.text(0.5, 0.95, f'Set Overlap: {overall_stats["set_overlap_evaluable"]:,} evaluable items\n'
                        f'Majority ID: {overall_stats["majority_id_evaluable"]:,} evaluable items',
             transform=ax2.transAxes, ha='center', va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
             fontsize=10)

    plt.tight_layout()

    # Save the plot
    plt.savefig(output_viz_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Visualization saved: {output_viz_path}")

    # Also save as PDF
    output_file_pdf = output_viz_path.with_suffix('.pdf')
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Visualization saved (PDF): {output_file_pdf}")

    plt.close()


def print_results(stats_df, breakdown_df, overall_stats):
    """Print formatted results"""
    print("\n" + "="*80)
    print("OVERALL MAPPING ACCURACY")
    print("="*80)

    print(f"Total Items: {overall_stats['total_items']:,}")
    print(f"Set Overlap Accuracy: {overall_stats['set_overlap_accurate']:,} / "
          f"{overall_stats['set_overlap_evaluable']:,} evaluable = "
          f"{overall_stats['set_overlap_accuracy_rate']:.1f}%")
    print(f"Majority ID Accuracy: {overall_stats['majority_id_accurate']:,} / "
          f"{overall_stats['majority_id_evaluable']:,} evaluable = "
          f"{overall_stats['majority_id_accuracy_rate']:.1f}%")

    print("\n" + "="*80)
    print("MAPPING ACCURACY BY CONFIDENCE CATEGORY")
    print("="*80)

    print(f"\n{'Confidence':<10} {'Total':<8} {'Set Overlap':<25} {'Majority ID':<25}")
    print(f"{'Category':<10} {'Items':<8} {'Evaluable':<10} {'Accurate':<10} {'Rate':<5} {'Evaluable':<10} {'Accurate':<10} {'Rate':<5}")
    print("-" * 80)

    for _, row in stats_df.iterrows():
        print(f"{row['confidence']:<10} {row['total_items']:<8} "
              f"{row['set_overlap_evaluable']:<10} {row['set_overlap_accurate']:<10} "
              f"{row['set_overlap_accuracy_rate']:<5.1f}% "
              f"{row['majority_id_evaluable']:<10} {row['majority_id_accurate']:<10} "
              f"{row['majority_id_accuracy_rate']:<5.1f}%")

    print("\n" + "="*60)
    print("DETAILED ACCURACY PATTERNS")
    print("="*60)

    print(f"\n{'Confidence':<10} {'Both':<8} {'Set Only':<10} {'Majority Only':<13} {'Neither':<10}")
    print(f"{'Category':<10} {'Accurate':<8} {'Accurate':<10} {'Accurate':<13} {'Accurate':<10}")
    print("-" * 60)

    for _, row in breakdown_df.iterrows():
        print(f"{row['confidence']:<10} {row['both_accurate']:<8} "
              f"{row['set_overlap_accurate_only']:<10} {row['majority_id_accurate_only']:<13} "
              f"{row['neither_accurate']:<10}")


def main():
    parser = argparse.ArgumentParser(description='Calculate mapping accuracy by confidence category')
    parser.add_argument('--input',
                       default='src/mapping/results/kraken/arivale_metabolites_v4_a_full_results.tsv',
                       help='Input TSV file path')
    parser.add_argument('--output',
                       default='src/mapping/results/kraken/',
                       help='Output directory for results')
    parser.add_argument('--no-viz', action='store_true',
                       help='Skip generating visualization')
    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.output)

    output_summary_path = output_dir / f"a_validation_{input_file.name.removesuffix('.tsv')}_summary.txt"
    output_viz_path = output_dir / f"a_validation_{input_file.name.removesuffix('.tsv')}_chart.png"

    if not input_file.exists():
        print(f"Error: Input file {input_file} not found!")
        return

    # Load data
    print(f"Loading data from: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    print(f"Loaded {len(df)} rows")
    print(f"Confidence categories: {df['confidence'].value_counts().to_dict()}")

    # Analyze accuracy
    stats_df, overall_stats = analyze_accuracy_by_confidence(df)
    breakdown_df = generate_detailed_breakdown(df)

    # Print results to console
    print_results(stats_df, breakdown_df, overall_stats)

    # Create visualization
    if not args.no_viz:
        print("\nCreating visualization...")
        create_accuracy_visualization(stats_df, overall_stats, output_viz_path)

    # Save results to file
    print(f"\nSaving results to: {output_summary_path}")
    with open(output_summary_path, 'w') as f:
        f.write("MAPPING ACCURACY BY CONFIDENCE CATEGORY\n")
        f.write("="*80 + "\n\n")

        f.write("OVERALL ACCURACY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total Items: {overall_stats['total_items']:,}\n")
        f.write(f"Set Overlap Accuracy: {overall_stats['set_overlap_accurate']:,} / "
               f"{overall_stats['set_overlap_evaluable']:,} evaluable = "
               f"{overall_stats['set_overlap_accuracy_rate']:.1f}%\n")
        f.write(f"Majority ID Accuracy: {overall_stats['majority_id_accurate']:,} / "
               f"{overall_stats['majority_id_evaluable']:,} evaluable = "
               f"{overall_stats['majority_id_accuracy_rate']:.1f}%\n\n")

        f.write("SUMMARY STATISTICS BY CONFIDENCE\n")
        f.write("-" * 40 + "\n")
        f.write(stats_df.to_string(index=False))
        f.write("\n\n")

        f.write("DETAILED ACCURACY PATTERNS\n")
        f.write("-" * 40 + "\n")
        f.write(breakdown_df.to_string(index=False))
        f.write("\n\n")

        f.write("ACCURACY DEFINITIONS:\n")
        f.write("- Set Overlap Accuracy: Any overlap between kg_canonical_ids_provided and kg_canonical_ids_assigned\n")
        f.write("- Majority ID Accuracy: Exact match between kg_majority_canonical_id_provided and kg_majority_canonical_id_assigned\n")
        f.write("- Evaluable: Both provided and assigned values are non-empty\n")

    print(f"Analysis complete!")


if __name__ == "__main__":
    main()