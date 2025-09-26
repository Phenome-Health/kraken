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
    'uniprot': 'uniprotkb',
    'gene_id': 'ensembl',
    'Labcorp LOINC ID': 'loinc',
    'Quest LOINC ID': 'loinc',
    'PubChem_CID': 'pubchem.compound',
    'ChEBI_ID': 'chebi',
    'HMDB_ID': 'hmdb',
    'LM_ID': 'lm',
    'KEGG_ID': 'kegg.compound',
    'INCHI_KEY': 'inchikey',
    'INCHIKEY': 'inchikey',
    'RefMet_ID': 'rm',
    'SMILES': 'smiles'
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
    # input_path = input_dir / 'metabolomics_metadata.tsv'
    input_path = input_dir / 'metabolites_ChemicalAnnotation-Table1.tsv'  # Latest from Lance, with refmet in it


    id_cols = ['CAS', 'KEGG', 'HMDB', 'PUBCHEM', 'INCHIKEY', 'SMILES'] + ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']

    df = pd.read_table(input_path, dtype={'PUBCHEM': str})
    df[id_cols] = df[id_cols].replace('-', np.nan)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_metabolites_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def extract_protein_ids(input_dir: Path, output_dir: Path):
    input_path = input_dir / 'proteomics_metadata.tsv'

    id_cols = ['uniprot', 'gene_id']
    other_cols = ['name', 'panel', 'gene_name', 'gene_description', 'gene_id', 'transcript_id', 'protein_id']

    df = pd.read_table(input_path, usecols=other_cols + id_cols, skiprows=13)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_proteins_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def extract_clinical_lab_ids(input_dir: Path, output_dir: Path):
    input_path = input_dir / 'chemistries_metadata.tsv'

    id_cols = ['Labcorp LOINC ID', 'Quest LOINC ID']
    other_cols = ['Name', 'Display Name', 'Labcorp ID', 'Labcorp Name', 'Labcorp LOINC Name', 'Quest ID', 'Quest Name']

    df = pd.read_table(input_path, usecols=other_cols + id_cols, skiprows=13)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_clinicallabs_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def extract_lipid_ids(input_dir: Path, output_dir: Path):
    input_path = input_dir / 'Arivale_lipidomics_metadata_REFMETANNOT.tsv'

    id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    other_cols = ['Input.name', 'Standardized.name', 'Formula', 'Exact.mass', 'Super.class', 'Main.class', 'Sub.class']

    df = pd.read_table(input_path, usecols=other_cols + id_cols)
    df[id_cols] = df[id_cols].replace('-', np.nan)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_lipids_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def extract_lipid_ids_norefmet(input_dir: Path, output_dir: Path):
    input_path = input_dir / 'lipids_ChemicalAnnotation-Table1.tsv'

    id_cols = ['HMDB', 'KEGG']
    other_cols = ['CHEM_ID', 'SUPER_PATHWAY', 'SUB_PATHWAY', 'PATHWAY_SORTORDER', 'CHEMICAL_NAME', 'PLOT_NAME', 'PLATFORM']

    df = pd.read_table(input_path, usecols=other_cols + id_cols)
    df[id_cols] = df[id_cols].replace('-', np.nan)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_lipids_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def extract_lipid_ids_all(input_dir: Path, output_dir: Path):
    input_path_source = input_dir / 'lipids_ChemicalAnnotation-Table1.tsv'
    input_path_refmet = input_dir / 'Arivale_lipidomics_metadata_REFMETANNOT.tsv'
    source_df = pd.read_table(input_path_source)
    refmet_df = pd.read_table(input_path_refmet)
    df = pd.merge(source_df, refmet_df, left_on='CHEMICAL_NAME', right_on='Input.name')

    id_cols = ['HMDB', 'KEGG'] + ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']

    df[id_cols] = df[id_cols].replace('-', np.nan)
    print(df)

    df[['curies', 'invalid_ids']] = df.apply(lambda row: get_curies(row, id_cols), axis=1, result_type='expand')
    print(df)

    output_path = output_dir / 'arivale_lipids_with_curies.tsv'
    df.to_csv(output_path, sep='\t', index=False)


def main():
    input_dir = PROJECT_ROOT_PATH / 'input_data' / 'mapping' / 'arivale'
    output_dir = PROJECT_ROOT_PATH / 'src' / 'mapping'

    extract_metabolite_ids(input_dir, output_dir)
    # extract_protein_ids(input_dir, output_dir)
    # extract_clinical_lab_ids(input_dir, output_dir)
    extract_lipid_ids_all(input_dir, output_dir)


if __name__ == "__main__":
    main()
