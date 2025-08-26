

import csv
import logging
from pathlib import Path

import jsonlines


def harmonize_umls(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    logging.info(f"Harmonizing UMLS: {input_file} -> {nodes_output}, {edges_output}")
    umls_infores = 'infores:umls'

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
                observable_node = {'id': loinc_observable_curie, 
                                    'equivalent_ids': [loinc_observable_curie], 
                                    'categories': ['biolink:ClinicalFinding'],
                                    'provided_by': [umls_infores]}
                nodes[observable_node['id']] = observable_node
            
            # Add a node representing the LOINC part captured in this row (e.g., compound measured), if one doesn't yet exist
            loinc_part_curie = f"LOINC:{loinc_part_id}"
            umls_part_curie = f"UMLS:{umls_part_id}"
            if umls_part_curie not in nodes:
                part_node = {'id': umls_part_curie,
                             'equivalent_ids': [umls_part_curie, loinc_part_curie],
                             'categories': ['biolink:NamedThing'],
                             'provided_by': [umls_infores]}
                nodes[part_node['id']] = part_node
            
            # Add an edge connecting the observable node to its parts
            edge = {'subject': loinc_observable_curie, 
                    'object': umls_part_curie, 
                    'predicate': 'biolink:related_to', 
                    'primary_knowledge_source': 'infores:umls'}
            edges.append(edge)
    
    logging.info(f"Saving {len(nodes)} nodes and {len(edges)} edges")
    with jsonlines.open(nodes_output, 'w') as writer:
        writer.write_all(list(nodes.values()))
    with jsonlines.open(edges_output, 'w') as writer:
        writer.write_all(edges)
