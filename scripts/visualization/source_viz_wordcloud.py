"""
A claude-generated script for visualizing the flow of sources to kraken (includes edge sources only).
"""
import json
import os
import sys
from collections import Counter
from pathlib import Path

from wordcloud import WordCloud
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from utils.logging_config import setup_logging

PROJECT_ROOT_PATH = Path(__file__).parents[1]


setup_logging()


def generate_wordcloud(sources, output_filepath,
                       width=1200, height=800, max_words=100):
    """
    Generate a word cloud from the list of sources.

    Args:
        sources (list): List of source names (already cleaned)
        output_filepath (str): Output image filepath
        width (int): Width of the word cloud image
        height (int): Height of the word cloud image
        max_words (int): Maximum number of words to display
    """

    # Count frequency of each source
    source_counts = Counter(sources)

    # Print some statistics
    print(f"Total sources: {len(sources)}")
    print(f"Unique sources: {len(source_counts)}")
    print("\nTop 10 most frequent sources:")
    for source, count in source_counts.most_common(10):
        print(f"  {source}: {count}")

    # Generate word cloud
    wordcloud = WordCloud(
        width=width,
        height=height,
        max_words=max_words,
        background_color='white',
        colormap='viridis',  # You can change this to other colormaps like 'plasma', 'cool', etc.
        relative_scaling=0.5,
        min_font_size=32,  # Minimum font size for readability
        max_font_size=42,  # Cap the maximum font size
        collocations=False,  # Allow words to be placed closer together
        prefer_horizontal=0.8,  # Prefer horizontal text (helps with packing)
        random_state=42
    ).generate_from_frequencies(source_counts)

    # Create the plot
    plt.figure(figsize=(width / 100, height / 100))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.tight_layout(pad=0)

    # Save the image
    plt.savefig(output_filepath, dpi=300, bbox_inches='tight')
    print(f"\nWord cloud saved as: {output_filepath}")

    # Show the plot
    plt.show()

    return source_counts


def main():
    with open(f"{PROJECT_ROOT_PATH}/scripts/primary_ks_edges.json", 'r') as edges_file:
        source_dag_edges = json.load(edges_file)

    sources = [source_edge[0] for source_edge in source_dag_edges]
    generate_wordcloud(sources, f"{PROJECT_ROOT_PATH}/scripts/wordcloud.png")


# Example usage:
if __name__ == "__main__":
    main()
