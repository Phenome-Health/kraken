#!/usr/bin/env python3
"""
Standalone utility for generating metagraphs from harmonized JSONL files
Can be run independently or as part of the pipeline
"""

import argparse
import sys
from pathlib import Path
import logging

from .metagraph import (
    generate_metagraph_for_source,
    generate_metagraph_streaming,
    save_metagraph,
    generate_metagraph_summary,
    create_cytoscape_metagraph,
    compare_metagraphs
)

def main():
    parser = argparse.ArgumentParser(description='Generate metagraph from harmonized JSONL files')
    parser.add_argument('nodes_file', type=Path, help='Path to nodes JSONL file')
    parser.add_argument('edges_file', type=Path, help='Path to edges JSONL file')
    parser.add_argument('-o', '--output-dir', type=Path, required=True, help='Output directory for metagraph')
    parser.add_argument('-s', '--source-name', help='Source name (default: inferred from path)')
    parser.add_argument('--summary', action='store_true', help='Print human-readable summary')
    parser.add_argument('--cytoscape', action='store_true', help='Generate Cytoscape visualization file')
    parser.add_argument('--min-edges', type=int, default=1, help='Minimum edge count for Cytoscape (default: 1)')
    parser.add_argument('--compare', nargs='+', type=Path, help='Compare with other metagraph files')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Validate inputs
    if not args.nodes_file.exists():
        logging.error(f"Nodes file not found: {args.nodes_file}")
        sys.exit(1)
    
    if not args.edges_file.exists():
        logging.error(f"Edges file not found: {args.edges_file}")
        sys.exit(1)
    
    # Generate metagraph
    source_name = args.source_name or args.nodes_file.parent.name
    
    try:
        # Generate main metagraph file
        metagraph_file = generate_metagraph_for_source(
            args.nodes_file, args.edges_file, args.output_dir, source_name
        )
        
        logging.info(f"Metagraph generated: {metagraph_file}")
        
        # Generate summary if requested
        if args.summary:
            stats = generate_metagraph_streaming(args.nodes_file, args.edges_file, source_name)
            summary = generate_metagraph_summary(stats)
            
            print("\n" + summary + "\n")
            
            # Save summary to file
            summary_file = args.output_dir / f"{source_name}_metagraph_summary.txt"
            with open(summary_file, 'w') as f:
                f.write(summary)
            logging.info(f"Summary saved: {summary_file}")
        
        # Generate Cytoscape file if requested
        if args.cytoscape:
            stats = generate_metagraph_streaming(args.nodes_file, args.edges_file, source_name)
            cytoscape_file = args.output_dir / f"{source_name}_metagraph_cytoscape.json"
            create_cytoscape_metagraph(stats, cytoscape_file, args.min_edges)
            logging.info(f"Cytoscape file generated: {cytoscape_file}")
        
        # Compare with other metagraphs if requested
        if args.compare:
            comparison_files = [metagraph_file] + list(args.compare)
            comparison_output = args.output_dir / f"{source_name}_metagraph_comparison.json"
            compare_metagraphs(comparison_files, comparison_output)
            logging.info(f"Comparison generated: {comparison_output}")
        
        logging.info("Metagraph generation complete!")
        
    except Exception as e:
        logging.error(f"Failed to generate metagraph: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()