"""
A claude-generated script for visualizing the flow of sources to kraken (includes edge sources only).
"""

import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))
from utils.constants import PRIMARY_KS
from utils.kg_io import stream_edges_from_jsonl
from utils.logging_config import setup_logging

PROJECT_ROOT_PATH = Path(__file__).parents[1]


setup_logging()


def create_dag_from_edges(edges):
    """Create a DAG from a list of edges."""
    G = nx.DiGraph()
    G.add_edges_from(edges)

    if not nx.is_directed_acyclic_graph(G):
        raise ValueError("The graph contains cycles and is not a DAG!")

    return G


def get_hierarchical_positions(G, horizontal_spacing=2.0, vertical_spacing=1.5):
    """Calculate positions for nodes in a hierarchical layout with better spacing."""
    topo_order = list(nx.topological_sort(G))

    # Calculate levels
    levels = {}
    for node in topo_order:
        if G.in_degree(node) == 0:
            levels[node] = 0
        else:
            levels[node] = max(levels[pred] for pred in G.predecessors(node)) + 1

    # Group nodes by level
    level_groups = defaultdict(list)
    for node, level in levels.items():
        level_groups[level].append(node)

    # Calculate positions with better spacing
    pos = {}
    max_level = max(levels.values())

    for level, nodes in level_groups.items():
        y = (max_level - level) * vertical_spacing

        # Better horizontal distribution
        if len(nodes) == 1:
            x_positions = [0]
        else:
            total_width = (len(nodes) - 1) * horizontal_spacing
            x_positions = np.linspace(-total_width / 2, total_width / 2, len(nodes))

        for i, node in enumerate(sorted(nodes)):
            pos[node] = (x_positions[i], y)

    return pos, levels


def visualize_dag_plotly(G, title="Interactive DAG Visualization", width=1000, height=700):
    """
    Create an interactive, stylistic DAG visualization using Plotly.
    """
    pos, levels = get_hierarchical_positions(G)

    # Create edge traces
    edge_x = []
    edge_y = []
    edge_info = []

    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_info.append(f"{edge[0]} → {edge[1]}")

    # Create node traces
    node_x = []
    node_y = []
    node_text = []
    node_colors = []
    node_sizes = []

    # Color scheme based on levels
    colors = px.colors.qualitative.Set3

    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(str(node))

        # Color by level
        level = levels[node]
        node_colors.append(colors[level % len(colors)])

        # Size by degree (influence)
        total_degree = G.in_degree(node) + G.out_degree(node)
        node_sizes.append(max(20, 15 + total_degree * 5))

    # Create the plot
    fig = go.Figure()

    # Add edges
    fig.add_trace(
        go.Scatter(
            x=edge_x,
            y=edge_y,
            line=dict(width=2, color="rgba(125, 125, 125, 0.8)"),
            hoverinfo="none",
            mode="lines",
            showlegend=False,
        )
    )

    # Add arrow heads for edges
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        # Calculate arrow position (80% along the edge)
        arrow_x = x0 + 0.8 * (x1 - x0)
        arrow_y = y0 + 0.8 * (y1 - y0)

        # Calculate arrow direction
        dx = x1 - x0
        dy = y1 - y0
        length = np.sqrt(dx**2 + dy**2)
        if length > 0:
            dx, dy = dx / length, dy / length

        fig.add_annotation(
            x=arrow_x,
            y=arrow_y,
            ax=arrow_x - dx * 0.1,
            ay=arrow_y - dy * 0.1,
            xref="x",
            yref="y",
            axref="x",
            ayref="y",
            arrowhead=2,
            arrowsize=1.5,
            arrowwidth=2,
            arrowcolor="rgba(125, 125, 125, 0.8)",
            showlegend=False,
        )

    # Add nodes
    fig.add_trace(
        go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            marker=dict(size=node_sizes, color=node_colors, line=dict(width=2, color="white"), symbol="circle"),
            text=node_text,
            textposition="middle center",
            textfont=dict(size=12, color="black", family="Arial Black"),
            hovertemplate="<b>%{text}</b><br>Level: %{customdata}<extra></extra>",
            customdata=[levels[node] for node in G.nodes()],
            showlegend=False,
        )
    )

    fig.update_layout(
        title=dict(text=title, x=0.5, font=dict(size=20, family="Arial Black")),
        showlegend=False,
        hovermode="closest",
        margin=dict(b=20, l=5, r=5, t=40),
        annotations=[
            dict(
                text="Drag nodes to rearrange • Hover for details",
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0.005,
                y=-0.002,
                xanchor="left",
                yanchor="bottom",
                font=dict(color="gray", size=12),
            )
        ],
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor="white",
        width=width,
        height=height,
    )

    return fig


def visualize_dag_matplotlib_stylish(G, title="Stylistic DAG Visualization", figsize=(14, 10), style="modern"):
    """
    Create a highly stylistic DAG visualization using matplotlib.
    """
    # Set style
    plt.style.use("default")
    if style == "dark":
        plt.style.use("dark_background")

    pos, levels = get_hierarchical_positions(G, horizontal_spacing=3, vertical_spacing=2)

    fig, ax = plt.subplots(figsize=figsize, facecolor="white" if style != "dark" else "black")

    # Color schemes
    if style == "modern":
        node_colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD", "#98D8C8"]
        edge_color = "#2C3E50"
        text_color = "white"
    elif style == "dark":
        node_colors = ["#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7", "#DDA0DD", "#98D8C8"]
        edge_color = "#BDC3C7"
        text_color = "white"
    else:  # pastel
        node_colors = ["#FFB3BA", "#BAFFC9", "#BAE1FF", "#FFFFBA", "#FFDFBA", "#E0BBE4", "#D4F0D4"]
        edge_color = "#34495E"
        text_color = "black"

    # Draw edges with style
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]

        # Draw edge with gradient effect
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="->",
                lw=2.5,
                color=edge_color,
                alpha=0.7,
                shrinkA=25,
                shrinkB=25,
                connectionstyle="arc3,rad=0.1",
            ),
        )

    # Draw nodes with fancy styling
    for node in G.nodes():
        x, y = pos[node]
        level = levels[node]

        # Node size based on connections
        total_degree = G.in_degree(node) + G.out_degree(node)
        size = 800 + total_degree * 200

        # Color by level
        color = node_colors[level % len(node_colors)]

        # Create fancy node
        if style == "modern":
            # Gradient-like effect with multiple circles
            for i, alpha in enumerate([0.3, 0.5, 0.7, 1.0]):
                circle_size = size * (1.2 - i * 0.05)
                circle = plt.Circle(
                    (x, y), np.sqrt(circle_size / np.pi) / 10, color=color, alpha=alpha, zorder=2 - i * 0.1
                )
                ax.add_patch(circle)
        else:
            # Standard styled circle
            circle = plt.Circle((x, y), np.sqrt(size / np.pi) / 10, color=color, alpha=0.9, zorder=2)
            ax.add_patch(circle)

            # Add border
            border = plt.Circle(
                (x, y), np.sqrt(size / np.pi) / 10, fill=False, edgecolor="white", linewidth=3, zorder=3
            )
            ax.add_patch(border)

        # Add text with better styling
        ax.text(
            x,
            y,
            str(node),
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color=text_color,
            zorder=4,
            bbox=(
                dict(boxstyle="round,pad=0.1", facecolor=color, alpha=0.8, edgecolor="none")
                if style != "modern"
                else None
            ),
        )

    # Styling
    ax.set_xlim(min(x for x, y in pos.values()) - 1, max(x for x, y in pos.values()) + 1)
    ax.set_ylim(min(y for x, y in pos.values()) - 1, max(y for x, y in pos.values()) + 1)
    ax.set_aspect("equal")
    ax.axis("off")

    # Title with style
    plt.suptitle(title, fontsize=20, fontweight="bold", color=text_color if style == "dark" else "black", y=0.95)

    # Add level indicators
    level_y_positions = {}
    for node, level in levels.items():
        x, y = pos[node]
        if level not in level_y_positions:
            level_y_positions[level] = y

    for level, y_pos in level_y_positions.items():
        ax.text(
            min(x for x, y in pos.values()) - 0.8,
            y_pos,
            f"Level {level}",
            fontsize=10,
            alpha=0.7,
            color=text_color if style == "dark" else "gray",
            verticalalignment="center",
        )

    plt.tight_layout()
    return fig


# Example usage
if __name__ == "__main__":

    root_edges = [
        ("rtx-kg2", "kraken"),
        ("spoke", "kraken"),
        ("umls", "kraken"),
        ("lipidmaps", "kraken"),
        ("refmet", "kraken"),
    ]
    # Load KG2 edge sources
    logging.info("Loading KG2 sources..")
    kg2_edges_path = f"{PROJECT_ROOT_PATH}/artifacts/harmonized/kg2/edges.jsonl"
    kg2_sources = set()
    for edge in stream_edges_from_jsonl(kg2_edges_path):
        primary_ks = edge[PRIMARY_KS]
        kg2_sources.add(primary_ks.removeprefix("infores:"))
    kg2_edges = [(source, "rtx-kg2") for source in kg2_sources]

    # Load SPOKE edge sources
    logging.info("Loading SPOKE sources..")
    spoke_edges_path = f"{PROJECT_ROOT_PATH}/artifacts/harmonized/spoke/edges.jsonl"
    spoke_sources = set()
    for edge in stream_edges_from_jsonl(spoke_edges_path):
        primary_ks = edge[PRIMARY_KS]
        spoke_sources.add(primary_ks.removeprefix("infores:"))
    spoke_edges = [(source, "spoke") for source in spoke_sources]

    # Create the DAG
    logging.info("Creating DAG...")
    edges = root_edges + kg2_edges + spoke_edges
    dag = create_dag_from_edges(edges)

    # Option 1: Interactive Plotly visualization
    logging.info("Creating interactive Plotly visualization...")
    fig_plotly = visualize_dag_plotly(dag, "ML Pipeline DAG - Interactive")
    fig_plotly.show()

    # Option 2: Stylistic matplotlib versions for presentations
    logging.info("Creating presentation-ready matplotlib visualizations...")

    # Modern style - vibrant and professional
    fig1 = visualize_dag_matplotlib_stylish(dag, "ML Pipeline DAG - Modern Style", style="modern")
    plt.show()

    # Professional style - corporate-friendly
    fig2 = visualize_dag_matplotlib_stylish(dag, "ML Pipeline DAG - Professional Style", style="professional")
    plt.show()

    # Pastel style - soft and elegant
    fig3 = visualize_dag_matplotlib_stylish(dag, "ML Pipeline DAG - Pastel Style", style="pastel")
    plt.show()

    # Save figures with transparent backgrounds for slides
    fig_plotly.write_html("dag_interactive.html")
    fig1.savefig("dag_modern.png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    fig2.savefig("dag_professional.png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    fig3.savefig("dag_pastel.png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")

    # Also save with transparent backgrounds
    fig1.savefig("dag_modern_transparent.png", dpi=300, bbox_inches="tight", facecolor="none", edgecolor="none")
    fig2.savefig("dag_professional_transparent.png", dpi=300, bbox_inches="tight", facecolor="none", edgecolor="none")
    fig3.savefig("dag_pastel_transparent.png", dpi=300, bbox_inches="tight", facecolor="none", edgecolor="none")

    logging.info("Visualizations saved!")
