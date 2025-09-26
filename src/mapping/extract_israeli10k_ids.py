import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Dict

import pandas as pd
import numpy as np
from jsonlines import jsonlines

PROJECT_ROOT_PATH = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT_PATH))

from src.utils.identifiers import IdentifierNorm

VOCAB_MAP = {
    'PubChem_CID': 'pubchem.compound',
    'ChEBI_ID': 'chebi',
    'HMDB_ID': 'hmdb',
    'LM_ID': 'lm',
    'KEGG_ID': 'kegg.compound',
    'INCHI_KEY': 'inchikey',
    'RefMet_ID': 'rm'
}

IDNORM = IdentifierNorm(biolink_version="4.2.5")


def get_curies(row: pd.Series, id_cols: List[str]) -> Tuple[List[str], Dict[str, List[str]]]:
    # First load IDs of different types, handling multiple IDs in one cell as necessary
    local_ids_map = {id_col: [local_id for local_id in re.split(';', row[id_col])]
                     for id_col in id_cols if pd.notnull(row[id_col])}

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

    return list(all_curies), dict(invalid_ids)


def extract_lipid_ids(input_dir: Path, output_dir: Path):
    # source_path = input_dir / 'lipids_ChemicalAnnotation-Table1.tsv'
    refmet_path = input_dir / 'Israeli10k_website_lipidomics_metadata_REFMETANNOT.tsv'
    # refmet_path = input_dir / 'israeli10k_lipidomics_metadata_REFMETANNOT.tsv'
    # source_df = pd.read_table(source_path)
    # refmet_df = pd.read_table(refmet_path)
    # df = pd.merge(source_df, refmet_df, left_on='Description', right_on='Input.name')

    df = pd.read_table(refmet_path)

    id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']

    df[id_cols] = df[id_cols].replace('-', np.nan)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'israeli10k_lipids_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def main():
    input_dir = PROJECT_ROOT_PATH / 'input_data' / 'mapping' / 'israeli10k'
    output_dir = PROJECT_ROOT_PATH / 'src' / 'mapping'

    extract_lipid_ids(input_dir, output_dir)


if __name__ == "__main__":
    main()
