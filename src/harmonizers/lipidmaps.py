import logging
from pathlib import Path

from rdkit import Chem
from biomapper2.core.normalizer import Normalizer

from ..utils.constants import *
from ..utils.kg_io import save_to_jsonl
from ..utils.general import create_node


def harmonize_lipidmaps(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    logging.info(f"Harmonizing LIPID MAPS: {input_file} -> {nodes_output}, {edges_output}")
    normalizer = Normalizer(biolink_version=biolink_version)
    attribute_prop_names = ['CATEGORY', 'MAIN_CLASS', 'SUB_CLASS', 'CLASS_LEVEL4', 'INCHI']
    nodes = dict()

    # Create a supplier object to read the file
    supplier = Chem.SDMolSupplier(input_file, removeHs=False)
    for molecule in supplier:

        # Some molecules might fail to load, so it's good practice to check
        if molecule is None:
            logging.warning(f"Lipid Maps molecule failed to load")
            continue

        properties = molecule.GetPropsAsDict()

        # Transform the 'canonical' ID for this node into standard curie form
        lm_curie_dict, _ = normalizer.get_curies({'LM_ID': properties['LM_ID']}, stop_on_invalid_id=True)
        lm_curie, lm_iri = next(iter(lm_curie_dict.items()))

        # Grab all xrefs and transform into standardized curies
        equivalent_ids = {lm_curie}
        equiv_id_properties = ['INCHI_KEY', 'PUBCHEM_CID', 'CHEBI_ID', 'KEGG_ID', 'HMDB_ID', 'SWISSLIPIDS_ID',
                               'LIPIDBANK_ID', 'PLANTFA_ID', 'SMILES']
        equiv_curies_dict, _ = normalizer.get_curies({prop_name: properties[prop_name]
                                                      for prop_name in equiv_id_properties if properties.get(prop_name)},
                                                     stop_on_invalid_id=False)
        if equiv_curies_dict:
            equivalent_ids = equivalent_ids | set(equiv_curies_dict)

        # Grab all names/synonyms
        name = properties.get('NAME')
        lm_synonyms = properties.get('SYNONYMS', '').split(';')
        other_synonyms = [name, properties.get('SYSTEMATIC_NAME'), properties.get('ABBREVIATION')]
        synonyms = {synonym.strip() for synonym in (lm_synonyms + other_synonyms) if synonym}

        # Put together our node
        node = create_node(
            curie=lm_curie,
            categories=['biolink:SmallMolecule'],
            equivalent_ids=list(equivalent_ids),
            provided_by=[LIPIDMAPS_CURIE],
            name=name,
            iri=lm_iri,
            synonyms=list(synonyms),
            chemical_formula=properties['FORMULA'],
            exact_mass=properties['EXACT_MASS'],
            attributes={LIPIDMAPS_CURIE: {attr_name: properties[attr_name]
                                          for attr_name in attribute_prop_names if attr_name in properties}}
        )

        nodes[node[ID]] = node

    logging.info(f"Saving {len(nodes)} lipid maps nodes")
    save_to_jsonl(nodes.values(), nodes_output, mode='w')
    save_to_jsonl([], edges_output, mode='w')  # Just create an empty edges file..

