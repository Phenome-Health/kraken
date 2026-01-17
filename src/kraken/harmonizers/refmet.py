# refmet.py
import logging
from pathlib import Path
from typing import Any, Dict

from biomapper2.core.normalizer import Normalizer

from kraken.utils.constants import REFMET_INFORES
from kraken.utils.kg_io import save_to_jsonl, load_csv_to_dict_list
from kraken.utils.biolink_client import BiolinkClient
from .base import BaseHarmonizer


class RefMetHarmonizer:
    """Harmonizer for RefMet CSV files - doesn't use base class due to unique format"""

    source_name = "refmet"
    source_infores = REFMET_INFORES

    attribute_props = {'super_class', 'main_class', 'sub_class'}
    equiv_id_props = {'pubchem_cid', 'chebi_id', 'hmdb_id', 'lipidmaps_id', 'kegg_id', 'inchi_key'}

    def __init__(self, biolink_client: BiolinkClient):
        self.biolink = biolink_client
        self.normalizer = Normalizer(biolink_version=biolink_client.version)

    def harmonize(
            self,
            input_file: Path,
            nodes_output: Path,
            edges_output: Path,
    ):
        logging.info(f"Harmonizing {self.source_name}: {input_file} -> {nodes_output}, {edges_output}")

        nodes = {}

        for row in load_csv_to_dict_list(input_file):
            node = self._harmonize_row(row)
            if node:
                nodes[node['id']] = node

        logging.info(f"Saving {len(nodes)} RefMet nodes")
        save_to_jsonl(nodes.values(), nodes_output, mode='w')
        save_to_jsonl([], edges_output, mode='w')  # Empty edges file

        logging.info(f"{self.source_name} harmonization complete: {len(nodes)} nodes, 0 edges")

    def _harmonize_row(self, row: Dict[str, Any]) -> Dict[str, Any] | None:
        # Transform the 'canonical' ID into standard curie form
        # Note: original has ' refmet_id' with leading space
        rm_curie_dict, _ = self.normalizer.get_curies(
            {'refmet': row[' refmet_id']},
            stop_on_invalid_id=True
        )
        rm_curie, rm_iri = next(iter(rm_curie_dict.items()))

        # Grab all xrefs and transform into standardized curies
        equivalent_ids = {rm_curie}
        equiv_curies_dict, _ = self.normalizer.get_curies(
            {prop: row[prop] for prop in self.equiv_id_props if row.get(prop)},
            stop_on_invalid_id=False
        )
        if equiv_curies_dict:
            equivalent_ids |= set(equiv_curies_dict)

        # Put together our node
        name = row['refmet_name']

        return BaseHarmonizer.create_node(
            source_infores=self.source_infores,
            curie=rm_curie,
            categories=['biolink:SmallMolecule'],
            equivalent_ids=list(equivalent_ids),
            provided_by=self.source_infores,
            name=name,
            urls=rm_iri,
            chemical_formula=row['formula'],
            exact_mass=row['exactmass'],
            attributes={k: row[k] for k in self.attribute_props}
        )
