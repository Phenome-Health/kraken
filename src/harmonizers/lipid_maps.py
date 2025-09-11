import logging
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors

from ..utils.constants import *
from ..utils.identifiers import IdentifierNorm


def harmonize_lipid_maps(input_file: Path, nodes_output: Path, edges_output: Path, biolink_version: str):
    logging.info(f"Harmonizing LIPID MAPS: {input_file} -> {nodes_output}, {edges_output}")

    id_norm = IdentifierNorm(biolink_version)

    # Create a supplier object to read the file
    supplier = Chem.SDMolSupplier(input_file, removeHs=False)
    for molecule in supplier:
        # Some molecules might fail to load, so it's good practice to check
        if molecule is None:
            logging.warning(f"Lipid Maps molecule failed to load")
            continue

        # Get all properties for the current molecule as a dictionary
        properties = molecule.GetPropsAsDict()



        lm_curie, lm_iri = id_norm.construct_curie(properties['LM_ID'], 'lm')
        print(properties['LM_ID'], lm_curie, lm_iri)
        if lm_curie:
            pass
        else:
            sys.exit(1)



    sys.exit(1)


