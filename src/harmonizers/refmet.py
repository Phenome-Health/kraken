import logging
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from ..utils.constants import *
from ..utils.identifiers import IdentifierNorm
from ..utils.kg_io import save_to_jsonl, load_csv_to_dict_list
from ..utils.general import create_node


def harmonize_refmet(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    logging.info(f"Harmonizing REFMET: {input_file} -> {nodes_output}, {edges_output}")
    id_norm = IdentifierNorm(biolink_version)
    nodes = dict()

    for row in load_csv_to_dict_list(input_file):

        rm_curie, rm_iri = id_norm.construct_curie(row[' refmet_id'], 'rm', stop_on_failure=True)
        # Grab all IDs/xrefs and transform into standardized curies
        equivalent_ids = {rm_curie}
        prefix_map = {
            'pubchem_cid': 'pubchem.compound',
            'chebi_id': 'chebi',
            'hmdb_id': 'hmdb',
            'lipidmaps_id': 'lm',
            'kegg_id': ['kegg.compound', 'kegg.drug'], # They give both of these in this field
            'inchi_key': 'inchikey'
        }
        for prop_name, prefix_entry in prefix_map.items():
            if row.get(prop_name):
                equiv_id = str(row[prop_name]).strip()
                if equiv_id:  # Sometimes stripping makes it an empty string (was only a space)
                    equiv_curie, _ = id_norm.construct_curie(equiv_id, prefix_entry, stop_on_failure=True)
                    if equiv_curie and equiv_curie != KNOWN_INVALID:
                        equivalent_ids.add(equiv_curie)

        # Put together our node
        name = row['refmet_name']
        node = create_node(curie=rm_curie,
                           categories=['biolink:SmallMolecule'],
                           equivalent_ids=list(equivalent_ids),
                           provided_by=[REFMET_CURIE],
                           name=name,
                           iri=rm_iri,
                           synonyms=[name],
                           chemical_formula=row['formula'],
                           exact_mass=row['exactmass'])

        # Tack on other attributes
        other_prop_names = ['super_class', 'main_class', 'sub_class']
        other_props = {other_prop_name: row[other_prop_name] for other_prop_name in other_prop_names}
        node['refmet_info'] = other_props

        nodes[node[ID]] = node

    logging.info(f"Saving {len(nodes)} refmet nodes")
    save_to_jsonl(nodes.values(), nodes_output, mode='w')
    save_to_jsonl([], edges_output, mode='w')  # Just create an empty edges file..

