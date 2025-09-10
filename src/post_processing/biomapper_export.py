"""
Biomapper export utilities
Exports node type-specific files for biomapper module
"""

import csv
from pathlib import Path
import json
import logging
from typing import Any, List
from collections import defaultdict
from ..utils.kg_io import stream_nodes_from_jsonl

FILE_MAP = {
    "biolink:BiologicalProcess": f"biological_processes",
    "biolink:CellularComponent": f"cellular_components",
    "biolink:ChemicalEntity": f"chemicals",
    "biolink:Disease": f"diseases",
    "biolink:Drug": f"drugs",
    "biolink:Gene": f"genes",
    "biolink:SmallMolecule": f"metabolites",  # TODO: Will this capture all? Some LOINC:LP nodes have type protein...
    "biolink:MolecularActivity": f"molecular_activities",
    "biolink:Pathway": f"pathways",
    "biolink:PhenotypicFeature": f"phenotypes",
    "biolink:Protein": f"proteins",
    "biolink:ClinicalFinding": f"clinical_findings"
}
ARRAY_DELIMITER = '||'


def write_to_csv(items: List[List[Any]], file_path: str, mode: str):
    with open(file_path, mode=mode) as tsv_file:
        writer = csv.writer(tsv_file)
        writer.writerows(items)


def export_for_biomapper(nodes_path: Path, output_dir: Path, kraken_version: str):
    file_paths = {category: f"{output_dir}/kraken_{kraken_version}_{file_core_name}.csv"
                  for category, file_core_name in FILE_MAP.items()}
    headers = ["id", "name", "category", "description", "synonyms", "xrefs"]
    for file_path in file_paths.values():
        write_to_csv([headers], file_path, 'w+')

    counter = 0
    node_counts_by_type = defaultdict(int)
    for node in stream_nodes_from_jsonl(nodes_path):
        for category in node['entity_types_ancestral']:
            if category in FILE_MAP:
                node_counts_by_type[category] += 1
                if any(ARRAY_DELIMITER in synonym for synonym in node.get('synonyms', [])):
                    raise ValueError(f"Found node with a pipe in one of its synonyms: {node}")
                synonyms_joined = ARRAY_DELIMITER.join(node.get('synonyms', []))
                equiv_ids_joined = ARRAY_DELIMITER.join(node.get('equivalent_ids', []))
                row = [node["id"], node.get("name"), category, node.get('description'), synonyms_joined, equiv_ids_joined]
                write_to_csv([row], file_paths[category], 'a')
    
    logging.info(f"Exported nodes in {len(node_counts_by_type)} categories")
    
    create_biomapper_summary(node_counts_by_type, output_dir)


def create_biomapper_summary(node_counts_by_type: dict, output_dir: Path):
    """Create a summary file with statistics about the export"""
    summary = {
        'export_summary': {
            'total_node_types': len(node_counts_by_type),
            'node_counts_by_type': node_counts_by_type,
            'total_nodes_exported': sum(node_counts_by_type.values())
        }
    }

    summary_path = output_dir / "export_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logging.info(f"Export summary saved to: {summary_path}")
