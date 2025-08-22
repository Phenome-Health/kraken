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
    "biolink:BiologicalProcess": "kg2c_biological_processes",
    "biolink:CellularComponent": "kg2c_cellular_components",
    "biolink:ChemicalEntity": "kg2c_chemicals",
    "biolink:Disease": "kg2c_diseases",
    "biolink:Drug": "kg2c_drugs",
    "biolink:Gene": "kg2c_genes",
    "biolink:SmallMolecule": "kg2c_metabolites",
    "biolink:MolecularActivity": "kg2c_molecular_activities",
    "biolink:Pathway": "kg2c_pathways",
    "biolink:PhenotypicFeature": "kg2c_phenotypes",
    "biolink:Protein": "kg2c_proteins"
}
ARRAY_DELIMITER = '||'


def write_to_csv(items: List[List[Any]], file_path: str, mode: str):
    with open(file_path, mode=mode) as tsv_file:
        writer = csv.writer(tsv_file)
        writer.writerows(items)


def export_for_biomapper(nodes_path: Path, output_dir: Path):
    headers = ["id", "name", "category", "description", "synonyms", "xrefs"]
    for file_name in FILE_MAP.values():
        write_to_csv([headers], f"{output_dir}/{file_name}.csv", 'w+')

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
                write_to_csv([row], f"{output_dir}/{FILE_MAP[category]}.csv", 'a')
        counter += 1
        if counter % 1000000 == 0:
            logging.info(f"Have processed {counter} nodes for biomapper export...")
    
    logging.info(f"Found nodes in {len(node_counts_by_type)} categories")
    
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
