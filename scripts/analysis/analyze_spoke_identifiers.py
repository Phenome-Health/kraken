#!/usr/bin/env python3
"""
Script to analyze SPOKE identifier patterns by gathering one example node
for each node type-source combination.

This helps understand how identifiers are formatted across different
node types and data sources in SPOKE.
"""

import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any

import jsonlines


def analyze_spoke_identifiers(input_file: Path, output_file: Path = None):
    """
    Analyze SPOKE data to gather example nodes for each type-source combination.

    Args:
        input_file: Path to SPOKE JSONL file
        output_file: Path to save results (optional, defaults to input_file_spoke_examples.json)
    """
    if output_file is None:
        output_file = input_file.parent / f"{input_file.stem}_spoke_examples.json"

    # Track combinations we've seen and store examples
    seen_combinations: set[tuple[str, str]] = set()
    examples: dict[str, dict[str, Any]] = defaultdict(dict)

    # Statistics
    total_nodes = 0

    logging.info(f"Analyzing SPOKE data from {input_file}")

    with jsonlines.open(input_file, "r") as reader:
        for line_num, record in enumerate(reader, 1):

            # Stop once we hit edges (nodes come first)
            if record.get("type") == "relationship":
                logging.info("Reached the first edge. Stopping.")
                break

            total_nodes += 1

            # Extract node type and source(s)
            labels = record.get("labels", [])
            properties = record.get("properties", {})

            # Handle both 'source' (string) and 'sources' (list)
            sources_list = []
            if "source" in properties and properties["source"]:
                sources_list.append(properties["source"])
            if "sources" in properties and properties["sources"]:
                if isinstance(properties["sources"], list):
                    sources_list.extend(properties["sources"])
                else:
                    sources_list.append(str(properties["sources"]))

            # Remove duplicates and use 'unknown' if no sources found
            sources_list = list(set(sources_list)) if sources_list else ["unknown"]

            # Handle multiple labels - use first one as primary type
            node_type = labels[0] if labels else "unknown"

            # Check each source for new combinations
            for source in sources_list:
                combination = (node_type, source)

                # If we haven't seen this combination, store it as example
                if combination not in seen_combinations:
                    seen_combinations.add(combination)

                    # Store example with relevant fields
                    example = {
                        "node_id": record.get("id"),
                        "labels": labels,
                        "identifier": properties.get("identifier"),
                        "name": properties.get("name"),
                        "all_sources": sources_list,
                        "property_names": sorted(list(set(properties.keys()).difference({"Mate_Version"}))),
                        "line_number": line_num,
                        "full_node": record,
                    }

                    examples[node_type][source] = example

                    logging.debug(f"Found new combination: {node_type} from {source}")

            if total_nodes % 100000 == 0:
                logging.info(f"Processed {total_nodes} nodes, found {len(seen_combinations)} combinations")

    # Create summary
    summary = {
        "total_nodes_processed": total_nodes,
        "unique_combinations": len(seen_combinations),
        "node_types": list(examples.keys()),
        "examples_by_type": {},
    }

    # Process examples for better organization
    for node_type, sources in examples.items():
        summary["examples_by_type"][node_type] = {
            "source_count": len(sources),
            "sources": list(sources.keys()),
            "examples": sources,
        }

    # Save results
    with open(output_file, "w") as f:
        json.dump({"summary": summary, "examples": examples}, f, indent=2)

    logging.info("Analysis complete!")
    logging.info(f"Found {len(seen_combinations)} unique type-source combinations")
    logging.info(f"Node types: {sorted(examples.keys())}")
    logging.info(f"Results saved to: {output_file}")

    return output_file, summary


def print_summary(examples_file: Path):
    """Print a concise summary of the identifier patterns found."""
    with open(examples_file) as f:
        data = json.load(f)

    summary = data["summary"]
    examples = data["examples"]

    print("\n=== SPOKE Identifier Pattern Analysis ===")
    print(f"Total combinations: {summary['unique_combinations']}")
    print(f"Node types: {len(summary['node_types'])}")

    print("\nNode Type Breakdown:")
    for node_type in sorted(examples.keys()):
        sources = examples[node_type]
        print(f"  {node_type}: {len(sources)} sources")
        for source, example in sources.items():
            identifier = example["identifier"] or "NULL"
            print(f"    {source}: {identifier}")

    print("\nIdentifier Pattern Examples:")
    for node_type in sorted(examples.keys()):
        print(f"\n{node_type}:")
        for source, example in examples[node_type].items():
            identifier = example["identifier"] or "NULL"
            name = example["name"][:40] + "..." if example["name"] and len(example["name"]) > 40 else example["name"]
            properties = (
                ", ".join(example.get("property_names", []))[:60] + "..."
                if len(", ".join(example.get("property_names", []))) > 60
                else ", ".join(example.get("property_names", []))
            )
            print(f"  {source:15} | {identifier:20} | {name:43} | {properties}")


def main():
    parser = argparse.ArgumentParser(description="Analyze SPOKE identifier patterns by type and source")
    parser.add_argument("input_file", type=Path, help="SPOKE JSONL file to analyze")
    parser.add_argument("-o", "--output", type=Path, help="Output JSON file (default: auto-generated)")
    parser.add_argument("--summary", action="store_true", help="Print summary to console after analysis")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Validate input
    if not args.input_file.exists():
        logging.error(f"Input file not found: {args.input_file}")
        return 1

    try:
        # Analyze the data
        output_file, summary_data = analyze_spoke_identifiers(args.input_file, args.output)

        # Print summary if requested
        if args.summary:
            print_summary(output_file)

        return 0

    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
