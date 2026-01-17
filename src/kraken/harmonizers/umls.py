# umls.py
import csv
import logging
from pathlib import Path
from typing import Any

from kraken.harmonizers.base import BaseHarmonizer
from kraken.utils.biolink_client import BiolinkClient
from kraken.utils.constants import ID, NOT_PROVIDED, ROOT_CATEGORY, ROOT_PREDICATE, UMLS_INFORES
from kraken.utils.kg_io import save_to_jsonl


class UMLSHarmonizer:
    """Harmonizer for UMLS TSV files - doesn't use base class due to unique format"""

    source_name = "umls"
    source_infores = UMLS_INFORES

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client

    def harmonize(
        self,
        input_file: Path,
        nodes_output: Path,
        edges_output: Path,
    ):
        logging.info(f"Harmonizing {self.source_name}: {input_file} -> {nodes_output}, {edges_output}")

        nodes = {}
        edges = []

        with open(input_file) as tsv_file:
            reader = csv.reader(tsv_file, delimiter="\t")
            next(reader)  # Skip the header row
            for row in reader:
                row_nodes, row_edge = self._harmonize_row(row, nodes)
                nodes.update(row_nodes)
                edges.append(row_edge)

        logging.info(f"Saving {len(nodes)} nodes and {len(edges)} edges")
        save_to_jsonl(nodes.values(), nodes_output, mode="w")
        save_to_jsonl(edges, edges_output, mode="w")

        logging.info(f"{self.source_name} harmonization complete: {len(nodes)} nodes, {len(edges)} edges")

    def _harmonize_row(self, row: list[str], existing_nodes: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
        loinc_observable_id, loinc_part_id, umls_part_id = row

        new_nodes = {}

        # Add a LOINC node for the observable (e.g., clinical lab test), if one doesn't yet exist
        loinc_observable_curie = f"LOINC:{loinc_observable_id}"
        if loinc_observable_curie not in existing_nodes:
            observable_node = BaseHarmonizer.create_node(
                source_infores=self.source_infores,
                curie=loinc_observable_curie,
                equivalent_ids=[loinc_observable_curie],
                categories=["biolink:ClinicalFinding"],
                provided_by=self.source_infores,
            )
            new_nodes[observable_node[ID]] = observable_node

        # Add a node representing the LOINC part captured in this row (e.g., compound measured), if one doesn't exist
        loinc_part_curie = f"LOINC:{loinc_part_id}"
        umls_part_curie = f"UMLS:{umls_part_id}"
        if umls_part_curie not in existing_nodes:
            part_node = BaseHarmonizer.create_node(
                source_infores=self.source_infores,
                curie=umls_part_curie,
                equivalent_ids=[umls_part_curie, loinc_part_curie],
                categories=[ROOT_CATEGORY],
                provided_by=self.source_infores,
            )
            new_nodes[part_node[ID]] = part_node

        # Add an edge connecting the observable node to its parts
        edge = BaseHarmonizer.create_edge(
            source_infores=self.source_infores,
            subject_id=loinc_observable_curie,
            object_id=umls_part_curie,
            predicate=ROOT_PREDICATE,
            primary_ks=self.source_infores,
            knowledge_level="knowledge_assertion",
            agent_type=NOT_PROVIDED,
        )

        return new_nodes, edge
