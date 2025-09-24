import os
import re
import sys
from collections import defaultdict
from email.policy import default
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import numpy as np
from jsonlines import jsonlines

PROJECT_ROOT_PATH = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT_PATH))

from src.utils.identifiers import IdentifierNorm

VOCAB_MAP = {
    'CAS': 'cas',
    'KEGG': ['kegg.compound', 'kegg.drug'],
    'HMDB': 'hmdb',
    'PUBCHEM': 'pubchem.compound'
}
ID_COLS = list(VOCAB_MAP.keys())

IDNORM = IdentifierNorm(biolink_version="4.2.5")

# KNOWN_INVALID_IDS = {'CAS': {'120K5305', '11/2/3483', '10-2005-9', '59-007'}}  # Record invalid IDs in the source data


def get_curies(row: pd.Series) -> Tuple[List[str], defaultdict[str, List[str]]]:
    # First load IDs of different types, handling multiple IDs in one cell as necessary
    local_ids_map = {id_col: [local_id for local_id in re.split('[,;]', row[id_col])]
                     for id_col in ID_COLS if pd.notnull(row[id_col])}

    all_curies = set()
    invalid_ids = defaultdict(list)
    for id_col, vocab_local_ids in local_ids_map.items():
        for local_id in vocab_local_ids:
            curie = IDNORM.construct_curie(local_id, VOCAB_MAP[id_col], stop_on_failure=False)[0]
            if curie:
                all_curies.add(curie)
            else:
                # It failed validation; record this
                invalid_ids[id_col].append(local_id)

    return list(all_curies), invalid_ids


def extract_metabolite_ids(file_name: str):
    file_path = PROJECT_ROOT_PATH / 'input_data' / 'mapping' / 'arivale' / file_name

    relevant_cols = ['CHEMICAL_ID', 'SUB_PATHWAY', 'SUPER_PATHWAY', 'BIOCHEMICAL_NAME']

    metabolite_col_name = 'BIOCHEMICAL_NAME'
    df = pd.read_table(file_path, dtype={'PUBCHEM': str}, usecols=relevant_cols + ID_COLS, skiprows=13)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(get_curies, axis=1, result_type='expand')
    print(df)

    output_path = PROJECT_ROOT_PATH / 'src' / 'mapping' / 'arivale_metabolites_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def main():
    metabolites_file_name = "metabolomics_metadata.tsv"

    extract_metabolite_ids(metabolites_file_name)


if __name__ == "__main__":
    main()
