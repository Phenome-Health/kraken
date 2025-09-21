import logging
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem

from ..utils.constants import *
from ..utils.identifiers import IdentifierNorm
from ..utils.kg_io import save_to_jsonl
from ..utils.general import create_node, get_canonical_smiles


def harmonize_lipidmaps(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    logging.info(f"Harmonizing LIPID MAPS: {input_file} -> {nodes_output}, {edges_output}")
    id_norm = IdentifierNorm(biolink_version)

    nodes = dict()

    # Create a supplier object to read the file
    supplier = Chem.SDMolSupplier(input_file, removeHs=False)
    for molecule in supplier:

        # Some molecules might fail to load, so it's good practice to check
        if molecule is None:
            logging.warning(f"Lipid Maps molecule failed to load")
            continue

        properties = molecule.GetPropsAsDict()

        # Grab all IDs/xrefs and transform into standardized curies
        lm_curie, lm_iri = id_norm.construct_curie(properties['LM_ID'], 'lm', stop_on_failure=True)
        equivalent_ids = {lm_curie}
        prefix_map = {
            'INCHI_KEY': 'inchikey',
            'PUBCHEM_CID': 'pubchem.compound',
            'CHEBI_ID': 'chebi',
            'KEGG_ID': ['kegg.compound', 'kegg.drug'],  # They give both of these in this field
            'HMDB_ID': 'hmdb',
            'SWISSLIPIDS_ID': 'slm',
            "LIPIDBANK_ID": 'lipidbank',
            "PLANTFA_ID": 'plantfa',
            'SMILES': 'smiles',
        }
        for prop_name, prefix_entry in prefix_map.items():
            if prop_name in properties:
                equiv_id = str(properties[prop_name])

                # If this is a SMILES string, convert it to canonical format
                if prop_name == 'SMILES':
                    canonical_smiles = get_canonical_smiles(equiv_id)
                    equiv_id = canonical_smiles if canonical_smiles else equiv_id  # Just use given if standardization failed

                equiv_curie, _ = id_norm.construct_curie(equiv_id, prefix_entry, stop_on_failure=True)
                if equiv_curie and equiv_curie != KNOWN_INVALID:
                    equivalent_ids.add(equiv_curie)

        # Grab all names/synonyms
        name = properties.get('NAME')
        lm_synonyms = properties.get('SYNONYMS', '').split(';')
        other_synonyms = [name, properties.get('SYSTEMATIC_NAME'), properties.get('ABBREVIATION')]
        synonyms = {synonym.strip() for synonym in (lm_synonyms + other_synonyms) if synonym}

        # Put together our node
        node = create_node(curie=lm_curie,
                           categories=['biolink:SmallMolecule'],
                           equivalent_ids=list(equivalent_ids),
                           provided_by=[LIPIDMAPS_CURIE],
                           name=name,
                           iri=lm_iri,
                           synonyms=list(synonyms),
                           chemical_formula=properties['FORMULA'],
                           exact_mass=properties['EXACT_MASS'])

        # Tack on other attributes
        other_prop_names = ['CATEGORY', 'MAIN_CLASS', 'SUB_CLASS', 'CLASS_LEVEL4', 'INCHI']
        other_props = {other_prop_name: properties[other_prop_name]
                       for other_prop_name in other_prop_names if other_prop_name in properties}
        node['lipidmaps_info'] = other_props

        nodes[node[ID]] = node

    logging.info(f"Saving {len(nodes)} lipid maps nodes")
    save_to_jsonl(nodes.values(), nodes_output, mode='w')
    save_to_jsonl([], edges_output, mode='w')  # Just create an empty edges file..

