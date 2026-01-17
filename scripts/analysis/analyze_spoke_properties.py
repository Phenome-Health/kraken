import argparse
import json
from collections import defaultdict
from pathlib import Path

from kraken.utils.kg_io import stream_mixed_jsonl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file_path", type=Path, help="SPOKE JSONL file to analyze")
    args = parser.parse_args()

    node_property_counts = defaultdict(int)
    edge_property_counts = defaultdict(int)
    node_property_examples = dict()
    edge_property_examples = dict()

    for item in stream_mixed_jsonl(args.input_file_path):
        item_type = item.get("type")
        properties = item.get("properties")

        for property_name, value in properties.items():
            if item_type == "node":
                node_property_counts[property_name] += 1
                if property_name not in node_property_examples:
                    node_property_examples[property_name] = value
            elif item_type == "relationship":
                edge_property_counts[property_name] += 1
                if property_name not in edge_property_examples:
                    edge_property_examples[property_name] = value

    node_counts_sorted = dict(sorted(node_property_counts.items(), key=lambda x: x[1], reverse=True))
    edge_counts_sorted = dict(sorted(edge_property_counts.items(), key=lambda x: x[1], reverse=True))

    with open("scripts/spoke_properties.json", "w+") as output_file:
        json.dump(
            {
                "nodes": node_counts_sorted,
                "edges": edge_counts_sorted,
                "nodes_examples": node_property_examples,
                "edges_examples": edge_property_examples,
            },
            output_file,
            indent=2,
        )


if __name__ == "__main__":
    main()
