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
    'CAS': 'cas',
    'KEGG': ['kegg.compound', 'kegg.drug'],
    'HMDB': 'hmdb',
    'PUBCHEM': 'pubchem.compound',
    'uniprot': 'uniprotkb'
}

IDNORM = IdentifierNorm(biolink_version="4.2.5")


def get_curies(row: pd.Series, id_cols: List[str]) -> Tuple[List[str], Dict[str, List[str]]]:
    # First load IDs of different types, handling multiple IDs in one cell as necessary
    local_ids_map = {id_col: [local_id for local_id in re.split('[,;]', row[id_col])]
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


def extract_metabolite_ids(input_dir: Path, output_dir: Path):
    input_path = input_dir / 'metabolomics_metadata.tsv'

    id_cols = ['CAS', 'KEGG', 'HMDB', 'PUBCHEM']
    other_cols = ['CHEMICAL_ID', 'SUB_PATHWAY', 'SUPER_PATHWAY', 'BIOCHEMICAL_NAME']

    df = pd.read_table(input_path, dtype={'PUBCHEM': str}, usecols=other_cols + id_cols, skiprows=13)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_metabolites_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def extract_protein_ids(input_dir: Path, output_dir: Path):
    input_path = input_dir / 'proteomics_metadata.tsv'

    id_cols = ['uniprot']
    other_cols = ['name', 'panel', 'gene_name', 'gene_description', 'gene_id', 'transcript_id', 'protein_id']

    df = pd.read_table(input_path, dtype={'PUBCHEM': str}, usecols=other_cols + id_cols, skiprows=13)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_proteins_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def main():
    input_dir = PROJECT_ROOT_PATH / 'input_data' / 'mapping' / 'arivale'
    output_dir = PROJECT_ROOT_PATH / 'src' / 'mapping'

    extract_metabolite_ids(input_dir, output_dir)
    extract_protein_ids(input_dir, output_dir)


if __name__ == "__main__":
    main()
