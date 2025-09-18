#!/usr/bin/env python3
"""
Analyze node categories by vocabulary in the KRAKEN knowledge graph.
Creates multiple visualizations: heatmap, bubble grid, and network graph.
"""

import json
import pandas as pd
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx


def extract_vocabulary_from_curie(curie: str) -> str:
    """Extract vocabulary prefix from CURIE (e.g., 'CHEBI:12345' -> 'CHEBI')"""
    if ':' not in curie:
        return 'unknown'
    return curie.split(':', 1)[0]


def clean_category(category: str) -> str:
    """Remove 'biolink:' prefix from category names"""
    return category.replace('biolink:', '')


def count_nodes_by_vocab_category(nodes_file: Path) -> pd.DataFrame:
    """
    Count nodes by vocabulary and category using id_prefixes.
    Returns DataFrame with vocab, category, and count columns.
    """
    vocab_category_counts = defaultdict(Counter)
    
    print(f"Processing {nodes_file}...")
    with open(nodes_file, 'r') as f:
        for i, line in enumerate(f):
            if i % 100000 == 0:
                print(f"  Processed {i:,} nodes...")
            
            try:
                node = json.loads(line.strip())
                id_prefixes = node.get('id_prefixes', [])
                categories = node.get('entity_types', [])  # Use entity_types for ArangoDB export
                
                # Count each category for each vocabulary prefix
                for prefix in id_prefixes:
                    for category in categories:
                        clean_cat = clean_category(category)
                        vocab_category_counts[prefix][clean_cat] += 1
                    
            except (json.JSONDecodeError, KeyError) as e:
                print(f"Error processing line {i}: {e}")
                continue
    
    print(f"  Processed {i+1:,} total nodes")
    
    # Convert to DataFrame
    data = []
    for vocab, category_counts in vocab_category_counts.items():
        for category, count in category_counts.items():
            data.append({'vocabulary': vocab, 'category': category, 'count': count})
    
    return pd.DataFrame(data)


def filter_data_for_visualization(df: pd.DataFrame, min_count: int = 50) -> pd.DataFrame:
    """Filter out vocabularies and categories with fewer than min_count total nodes"""
    # Calculate totals by vocabulary and category
    vocab_totals = df.groupby('vocabulary')['count'].sum()
    category_totals = df.groupby('category')['count'].sum()
    
    # Filter vocabularies and categories with at least min_count
    valid_vocabs = vocab_totals[vocab_totals >= min_count].index
    valid_categories = category_totals[category_totals >= min_count].index
    
    # Filter the dataframe
    filtered_df = df[
        (df['vocabulary'].isin(valid_vocabs)) & 
        (df['category'].isin(valid_categories))
    ].copy()
    
    print(f"Filtered from {df['vocabulary'].nunique()} to {filtered_df['vocabulary'].nunique()} vocabularies")
    print(f"Filtered from {df['category'].nunique()} to {filtered_df['category'].nunique()} categories")
    
    return filtered_df


def create_heatmap(df: pd.DataFrame, output_dir: Path):
    """Create heatmap visualization"""
    # Filter data for better visualization
    df_filtered = filter_data_for_visualization(df, min_count=1000)
    
    # Pivot data for heatmap
    heatmap_data = df_filtered.pivot(index='vocabulary', columns='category', values='count').fillna(0)
    
    # Create heatmap
    plt.figure(figsize=(20, 12))
    sns.heatmap(heatmap_data, 
                annot=True, 
                fmt='g',
                cmap='YlOrRd',
                cbar_kws={'label': 'Node Count'},
                square=False)
    plt.title('Node Categories by Vocabulary - Heatmap (Filtered: ≥1000 counts)')
    plt.xlabel('Category')
    plt.ylabel('Vocabulary')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    output_file = output_dir / 'vocabulary_category_heatmap.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved heatmap: {output_file}")


def create_bubble_chart(df: pd.DataFrame, output_dir: Path):
    """Create bubble grid chart visualization"""
    # Filter data for better visualization
    df_filtered = filter_data_for_visualization(df, min_count=1000)
    
    # Add log-scaled count for color and size mapping
    import numpy as np
    df_viz = df_filtered.copy()
    df_viz['sqrt_count'] = df_viz['count'].apply(lambda x: max(1, x)).apply(lambda x: x**0.5)  # Square root scale for size
    df_viz['log_count'] = df_viz['count'].apply(lambda x: max(1, x)).apply(np.log10)  # Log10 scale for color
    
    # Create bubble chart with log scaling
    fig = px.scatter(df_viz, 
                     x='category', 
                     y='vocabulary',
                     size='sqrt_count',
                     color='log_count',  # Use log scale for color
                     hover_data=['count'],
                     title='Node Categories by Vocabulary - Bubble Chart (Filtered: ≥1000 counts)',
                     color_continuous_scale='Viridis',
                     size_max=50)
    
    # Calculate appropriate height based on number of vocabularies
    num_vocabularies = df_filtered['vocabulary'].nunique()
    height = max(1200, num_vocabularies * 25)  # At least 25 pixels per vocabulary
    
    fig.update_layout(
        width=1400,
        height=height,
        xaxis_tickangle=45,
        xaxis_title='Category',
        yaxis_title='Vocabulary',
        yaxis=dict(categoryorder='total ascending')  # Sort vocabularies by total count
    )
    
    # Update colorbar to show log scale
    fig.update_coloraxes(colorbar_title="Count (log₁₀ scale)")
    
    output_file = output_dir / 'vocabulary_category_bubble.html'
    fig.write_html(output_file)
    print(f"Saved bubble chart: {output_file}")


def create_network_graph(df: pd.DataFrame, output_dir: Path):
    """Create network graph visualization"""
    # Filter data for better visualization
    df_filtered = filter_data_for_visualization(df, min_count=1000)
    
    # Create network graph
    G = nx.Graph()
    
    # Add nodes
    vocabularies = df_filtered['vocabulary'].unique()
    categories = df_filtered['category'].unique()
    
    # Add vocabulary nodes (left side)
    for vocab in vocabularies:
        G.add_node(f"vocab_{vocab}", node_type='vocabulary', label=vocab, size=20)
    
    # Add category nodes (right side)  
    for cat in categories:
        G.add_node(f"cat_{cat}", node_type='category', label=cat, size=15)
    
    # Add edges with weights
    for _, row in df_filtered.iterrows():
        vocab_node = f"vocab_{row['vocabulary']}"
        cat_node = f"cat_{row['category']}"
        weight = row['count']
        G.add_edge(vocab_node, cat_node, weight=weight)
    
    # Create layout
    pos = {}
    vocab_nodes = [n for n in G.nodes() if n.startswith('vocab_')]
    cat_nodes = [n for n in G.nodes() if n.startswith('cat_')]
    
    # Position vocabulary nodes on left
    for i, node in enumerate(vocab_nodes):
        pos[node] = (-1, i - len(vocab_nodes)/2)
    
    # Position category nodes on right
    for i, node in enumerate(cat_nodes):
        pos[node] = (1, i - len(cat_nodes)/2)
    
    # Create interactive plot with plotly
    edge_x = []
    edge_y = []
    edge_weights = []
    
    for edge in G.edges(data=True):
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_weights.append(edge[2]['weight'])
    
    # Normalize edge weights for line width
    max_weight = max(edge_weights) if edge_weights else 1
    normalized_weights = [w/max_weight * 10 for w in edge_weights]
    
    fig = go.Figure()
    
    # Add edges
    for i in range(0, len(edge_x), 3):
        if i//3 < len(edge_weights):
            fig.add_trace(go.Scatter(
                x=edge_x[i:i+2],
                y=edge_y[i:i+2],
                mode='lines',
                line=dict(width=normalized_weights[i//3], color='rgba(125,125,125,0.5)'),
                hoverinfo='none',
                showlegend=False
            ))
    
    # Add vocabulary nodes
    vocab_x = [pos[node][0] for node in vocab_nodes]
    vocab_y = [pos[node][1] for node in vocab_nodes]
    vocab_labels = [G.nodes[node]['label'] for node in vocab_nodes]
    
    fig.add_trace(go.Scatter(
        x=vocab_x, y=vocab_y,
        mode='markers+text',
        marker=dict(size=20, color='lightblue'),
        text=vocab_labels,
        textposition='middle left',
        name='Vocabularies',
        hovertemplate='Vocabulary: %{text}<extra></extra>'
    ))
    
    # Add category nodes
    cat_x = [pos[node][0] for node in cat_nodes]
    cat_y = [pos[node][1] for node in cat_nodes]
    cat_labels = [G.nodes[node]['label'] for node in cat_nodes]
    
    fig.add_trace(go.Scatter(
        x=cat_x, y=cat_y,
        mode='markers+text',
        marker=dict(size=15, color='lightcoral'),
        text=cat_labels,
        textposition='middle right',
        name='Categories',
        hovertemplate='Category: %{text}<extra></extra>'
    ))
    
    fig.update_layout(
        title='Vocabulary-Category Network Graph (Filtered: ≥1000 counts)<br><sub>Edge thickness represents node count</sub>',
        showlegend=True,
        width=1400,
        height=1000,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='white'
    )
    
    output_file = output_dir / 'vocabulary_category_network.html'
    fig.write_html(output_file)
    print(f"Saved network graph: {output_file}")


def generate_summary_stats(df: pd.DataFrame, output_dir: Path):
    """Generate summary statistics"""
    print("\n=== SUMMARY STATISTICS ===")
    print(f"Total vocabularies: {df['vocabulary'].nunique()}")
    print(f"Total categories: {df['category'].nunique()}")
    print(f"Total node count: {df['count'].sum():,}")
    
    print(f"\nTop 10 vocabularies by node count:")
    vocab_totals = df.groupby('vocabulary')['count'].sum().sort_values(ascending=False)
    for vocab, count in vocab_totals.head(10).items():
        print(f"  {vocab}: {count:,}")
    
    print(f"\nTop 10 categories by node count:")
    cat_totals = df.groupby('category')['count'].sum().sort_values(ascending=False)
    for cat, count in cat_totals.head(10).items():
        print(f"  {cat}: {count:,}")
    
    # Save summary to file
    summary_file = output_dir / 'vocabulary_category_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("=== VOCABULARY-CATEGORY ANALYSIS SUMMARY ===\n\n")
        f.write(f"Total vocabularies: {df['vocabulary'].nunique()}\n")
        f.write(f"Total categories: {df['category'].nunique()}\n") 
        f.write(f"Total node count: {df['count'].sum():,}\n\n")
        
        f.write("Top vocabularies by node count:\n")
        for vocab, count in vocab_totals.items():
            f.write(f"  {vocab}: {count:,}\n")
            
        f.write(f"\nTop categories by node count:\n")
        for cat, count in cat_totals.items():
            f.write(f"  {cat}: {count:,}\n")
    
    print(f"Saved summary: {summary_file}")


def main():
    # File paths - use ArangoDB export
    nodes_file = Path("artifacts/export/arango/kraken_1.0.0_nodes_arango.jsonl")
    output_dir = Path("scripts")
    
    if not nodes_file.exists():
        print(f"Error: {nodes_file} not found!")
        return
    
    # Count nodes by vocabulary and category
    df = count_nodes_by_vocab_category(nodes_file)
    
    if df.empty:
        print("No data found!")
        return
    
    # Save raw data
    csv_file = output_dir / 'vocabulary_category_counts.csv'
    df.to_csv(csv_file, index=False)
    print(f"Saved raw data: {csv_file}")
    
    # Generate visualizations
    create_heatmap(df, output_dir)
    create_bubble_chart(df, output_dir)
    create_network_graph(df, output_dir)
    
    # Generate summary statistics
    generate_summary_stats(df, output_dir)


if __name__ == "__main__":
    main()