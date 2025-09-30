#!/usr/bin/env python3
"""
Visualize knowledge_level and agent_type counts from metagraph JSON files.
Creates bar charts showing the distribution of edge metadata attributes.
"""

import json
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from collections import Counter


def load_metagraph_data(json_file):
    """Load metagraph data from JSON file"""
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data


def create_metadata_visualization(data, output_dir, source_name):
    """Create visualization of knowledge_level and agent_type distributions"""

    # Extract metadata counts and total edges
    knowledge_levels = data.get('knowledge_levels', {})
    agent_types = data.get('agent_types', {})
    primary_knowledge_sources = data.get('primary_knowledge_sources', {})
    total_edges = data.get('total_edges', 0)

    # Determine layout based on available data (excluding primary knowledge sources)
    available_plots = []
    if knowledge_levels:
        available_plots.append('knowledge_levels')
    if agent_types:
        available_plots.append('agent_types')

    if not available_plots:
        print("No metadata found to visualize")
        return

    # Set up the plot - simple layout for knowledge levels and agent types only
    n_plots = len(available_plots)
    if n_plots == 1:
        fig, ax = plt.subplots(1, 1, figsize=(14, 10))
        axes = [ax]
    elif n_plots == 2:
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    else:
        # This shouldn't happen since we only have 2 plot types now
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))

    plot_idx = 0

    # Plot knowledge levels
    if 'knowledge_levels' in available_plots:
        ax = axes[plot_idx]
        plot_idx += 1

        # Sort by count descending
        sorted_levels = sorted(knowledge_levels.items(), key=lambda x: x[1], reverse=True)
        levels, counts = zip(*sorted_levels) if sorted_levels else ([], [])

        # Calculate percentages based on sum of this metadata category
        category_total = sum(counts) if counts else 1
        percentages = [(count / category_total * 100) if category_total > 0 else 0 for count in counts]

        bars = ax.bar(range(len(levels)), percentages, color='lightblue', alpha=0.8)

        # Add value labels on bars with percentages and counts
        for i, bar in enumerate(bars):
            height = bar.get_height()
            count = counts[i]
            ax.annotate(f'$\\mathbf{{{height:.1f}\\%}}$\n({count:,})',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=14)

        ax.set_xlabel('Knowledge Level', fontweight='bold', fontsize=16)
        ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=16)
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels(levels, rotation=45, ha='right', fontsize=14)
        ax.set_ylim(0, max(100, max(percentages) * 1.1) if percentages else 100)
        ax.grid(True, alpha=0.3, axis='y')

    # Plot agent types
    if 'agent_types' in available_plots:
        ax = axes[plot_idx]
        plot_idx += 1

        # Sort by count descending
        sorted_agents = sorted(agent_types.items(), key=lambda x: x[1], reverse=True)
        agents, counts = zip(*sorted_agents) if sorted_agents else ([], [])

        # Calculate percentages based on sum of this metadata category
        category_total = sum(counts) if counts else 1
        percentages = [(count / category_total * 100) if category_total > 0 else 0 for count in counts]

        bars = ax.bar(range(len(agents)), percentages, color='lightcoral', alpha=0.8)

        # Add value labels on bars with percentages and counts
        for i, bar in enumerate(bars):
            height = bar.get_height()
            count = counts[i]
            ax.annotate(f'$\\mathbf{{{height:.1f}\\%}}$\n({count:,})',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=14)

        ax.set_xlabel('Agent Type', fontweight='bold', fontsize=16)
        ax.set_ylabel('Percentage (%)', fontweight='bold', fontsize=16)
        ax.set_xticks(range(len(agents)))
        ax.set_xticklabels(agents, rotation=45, ha='right', fontsize=14)
        ax.set_ylim(0, max(100, max(percentages) * 1.1) if percentages else 100)
        ax.grid(True, alpha=0.3, axis='y')



    plt.tight_layout()
    plt.subplots_adjust(top=0.85, bottom=0.15)  # Standard spacing for single row

    # Save the plot
    output_file = output_dir / f'{source_name.lower()}_metagraph_metadata_distribution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Visualization saved: {output_file}")

    # Also save as PDF
    output_file_pdf = output_dir / f'{source_name.lower()}_metagraph_metadata_distribution.pdf'
    plt.savefig(output_file_pdf, bbox_inches='tight', facecolor='white')
    print(f"Visualization saved (PDF): {output_file_pdf}")

    plt.close()


def print_summary_stats(data, source_name):
    """Print summary statistics to console"""
    print(f"\n{'='*60}")
    print(f"METAGRAPH METADATA SUMMARY: {source_name}")
    print(f"{'='*60}")

    total_edges = data.get('total_edges', 0)
    print(f"Total Edges: {total_edges:,}")

    # Knowledge levels
    knowledge_levels = data.get('knowledge_levels', {})
    if knowledge_levels:
        print(f"\nKnowledge Levels ({len(knowledge_levels)} distinct):")
        sorted_levels = sorted(knowledge_levels.items(), key=lambda x: x[1], reverse=True)
        for level, count in sorted_levels:
            percentage = (count / total_edges * 100) if total_edges > 0 else 0
            print(f"  {level}: {count:,} ({percentage:.1f}%)")

    # Agent types
    agent_types = data.get('agent_types', {})
    if agent_types:
        print(f"\nAgent Types ({len(agent_types)} distinct):")
        sorted_agents = sorted(agent_types.items(), key=lambda x: x[1], reverse=True)
        for agent, count in sorted_agents:
            percentage = (count / total_edges * 100) if total_edges > 0 else 0
            print(f"  {agent}: {count:,} ({percentage:.1f}%)")

    # Primary knowledge sources
    primary_knowledge_sources = data.get('primary_knowledge_sources', {})
    if primary_knowledge_sources:
        print(f"\nPrimary Knowledge Sources ({len(primary_knowledge_sources)} distinct):")
        sorted_sources = sorted(primary_knowledge_sources.items(), key=lambda x: x[1], reverse=True)
        for source, count in sorted_sources:
            percentage = (count / total_edges * 100) if total_edges > 0 else 0
            print(f"  {source}: {count:,} ({percentage:.1f}%)")


def save_summary_report(data, source_name, output_dir):
    """Save summary report to text file"""
    output_file = output_dir / f'{source_name.lower()}_metagraph_metadata_summary.txt'

    with open(output_file, 'w') as f:
        f.write(f"METAGRAPH METADATA SUMMARY: {source_name}\n")
        f.write("="*60 + "\n\n")

        total_edges = data.get('total_edges', 0)
        f.write(f"Total Edges: {total_edges:,}\n")

        # Knowledge levels
        knowledge_levels = data.get('knowledge_levels', {})
        if knowledge_levels:
            f.write(f"\nKnowledge Levels ({len(knowledge_levels)} distinct):\n")
            sorted_levels = sorted(knowledge_levels.items(), key=lambda x: x[1], reverse=True)
            for level, count in sorted_levels:
                percentage = (count / total_edges * 100) if total_edges > 0 else 0
                f.write(f"  {level}: {count:,} ({percentage:.1f}%)\n")

        # Agent types
        agent_types = data.get('agent_types', {})
        if agent_types:
            f.write(f"\nAgent Types ({len(agent_types)} distinct):\n")
            sorted_agents = sorted(agent_types.items(), key=lambda x: x[1], reverse=True)
            for agent, count in sorted_agents:
                percentage = (count / total_edges * 100) if total_edges > 0 else 0
                f.write(f"  {agent}: {count:,} ({percentage:.1f}%)\n")

        # Primary knowledge sources
        primary_knowledge_sources = data.get('primary_knowledge_sources', {})
        if primary_knowledge_sources:
            f.write(f"\nPrimary Knowledge Sources ({len(primary_knowledge_sources)} distinct):\n")
            sorted_sources = sorted(primary_knowledge_sources.items(), key=lambda x: x[1], reverse=True)
            for source, count in sorted_sources:
                percentage = (count / total_edges * 100) if total_edges > 0 else 0
                f.write(f"  {source}: {count:,} ({percentage:.1f}%)\n")

    print(f"Summary report saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize metagraph edge metadata distributions')
    parser.add_argument('json_file', help='Path to metagraph JSON file')
    parser.add_argument('--output', default='artifacts/metagraphs',
                       help='Output directory (default: artifacts/metagraphs)')
    parser.add_argument('--no-viz', action='store_true',
                       help='Skip generating visualization')
    args = parser.parse_args()

    json_file = Path(args.json_file)
    output_dir = Path(args.output)

    if not json_file.exists():
        print(f"Error: JSON file {json_file} not found!")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading metagraph data from: {json_file}")
    data = load_metagraph_data(json_file)

    # Extract source name from filename or data
    source_name = data.get('source', json_file.stem.replace('_metagraph', ''))

    # Print summary statistics
    print_summary_stats(data, source_name)

    # Create visualization
    if not args.no_viz:
        print(f"\nCreating visualization...")
        create_metadata_visualization(data, output_dir, source_name)

    # Save summary report
    print(f"\nSaving summary report...")
    save_summary_report(data, source_name, output_dir)

    print(f"\nAnalysis complete!")

    print(f"\nUsage examples:")
    print(f"  python src/utils/visualize_metagraph_metadata.py artifacts/metagraphs/kraken_metagraph.json")
    print(f"  python src/utils/visualize_metagraph_metadata.py artifacts/metagraphs/kg2_metagraph.json --output custom_output/")


if __name__ == "__main__":
    main()