import logging
from pathlib import Path

from biomapper2.core.normalizer import Normalizer

from ..utils.constants import *
from ..utils.kg_io import save_to_jsonl, load_csv_to_dict_list
from ..utils.general import create_node


def harmonize_refmet(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    logging.info(f"Harmonizing REFMET: {input_file} -> {nodes_output}, {edges_output}")
    normalizer = Normalizer(biolink_version=biolink_version)
    attribute_prop_names = ['super_class', 'main_class', 'sub_class']
    nodes = dict()

    for row in load_csv_to_dict_list(input_file):

        # Transform the 'canonical' ID for this node into standard curie form
        rm_curie_dict, _ = normalizer.get_curies({'refmet': row[' refmet_id']}, stop_on_invalid_id=True)
        rm_curie, rm_iri = next(iter(rm_curie_dict.items()))

        # Grab all xrefs and transform into standardized curies
        equivalent_ids = {rm_curie}
        equiv_id_properties = ['pubchem_cid', 'chebi_id', 'hmdb_id', 'lipidmaps_id', 'kegg_id', 'inchi_key']
        equiv_curies_dict, _ = normalizer.get_curies({prop_name: row[prop_name]
                                                      for prop_name in equiv_id_properties if row.get(prop_name)},
                                                     stop_on_invalid_id=False)
        if equiv_curies_dict:
            equivalent_ids = equivalent_ids | set(equiv_curies_dict)

        # Put together our node
        name = row['refmet_name']
        node = create_node(
            curie=rm_curie,
            categories=['biolink:SmallMolecule'],
            equivalent_ids=list(equivalent_ids),
            provided_by=[REFMET_CURIE],
            name=name,
            iri=rm_iri,
            synonyms=[name],
            chemical_formula=row['formula'],
            exact_mass=row['exactmass'],
            attributes={REFMET_CURIE: {attr_name: row[attr_name] for attr_name in attribute_prop_names}}
        )

        nodes[node[ID]] = node

    logging.info(f"Saving {len(nodes)} refmet nodes")
    save_to_jsonl(nodes.values(), nodes_output, mode='w')
    save_to_jsonl([], edges_output, mode='w')  # Just create an empty edges file..

