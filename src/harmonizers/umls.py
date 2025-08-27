

import csv
import logging
from pathlib import Path

import jsonlines
from ..utils.constants import *
from ..utils.kg_io import save_to_jsonl


def harmonize_umls(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str, build_metagraph: bool):
    logging.info(f"Harmonizing UMLS: {input_file} -> {nodes_output}, {edges_output}")

    nodes = dict()
    edges = []

    with open(input_file, 'r') as tsv_file:
        reader = csv.reader(tsv_file, delimiter='\t')
        next(reader)  # Skip the header row
        for row in reader:
            loinc_observable_id, loinc_part_id, umls_part_id = row

            # Add a LOINC node for the observable (e.g., clinical lab test), if one doesn't yet exist
            loinc_observable_curie = f"LOINC:{loinc_observable_id}"
            if loinc_observable_curie not in nodes:
                observable_node = {ID: loinc_observable_curie,
                                   EQUIVALENT_IDS: [loinc_observable_curie],
                                   CATEGORIES: ['biolink:ClinicalFinding'],
                                   PROVIDED_BY: [UMLS_INFORES]}
                nodes[observable_node[ID]] = observable_node
            
            # Add a node representing the LOINC part captured in this row (e.g., compound measured), if one doesn't yet exist
            loinc_part_curie = f"LOINC:{loinc_part_id}"
            umls_part_curie = f"UMLS:{umls_part_id}"
            if umls_part_curie not in nodes:
                part_node = {ID: umls_part_curie,
                             EQUIVALENT_IDS: [umls_part_curie, loinc_part_curie],
                             CATEGORIES: [ROOT_CATEGORY],
                             PROVIDED_BY: [UMLS_INFORES]}
                nodes[part_node[ID]] = part_node
            
            # Add an edge connecting the observable node to its parts
            edge = {SUBJECT: loinc_observable_curie,
                    OBJECT: umls_part_curie,
                    PREDICATE: ROOT_PREDICATE,
                    PRIMARY_KS: UMLS_INFORES}
            edges.append(edge)
    
    logging.info(f"Saving {len(nodes)} nodes and {len(edges)} edges")
    save_to_jsonl(nodes.values(), nodes_output, mode='w')
    save_to_jsonl(edges, edges_output, mode='w')
