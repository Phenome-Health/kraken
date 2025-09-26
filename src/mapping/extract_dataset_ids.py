import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import pandas as pd
import numpy as np
from jsonlines import jsonlines

PROJECT_ROOT_PATH = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT_PATH))

from src.utils.identifiers import IdentifierNorm

VOCAB_MAP = {
    'CAS': 'cas',
    'ChEBI_ID': 'chebi',
    'gene_id': 'ensembl',
    'HMDB': 'hmdb',
    'HMDB_ID': 'hmdb',
    'INCHI_KEY': 'inchikey',
    'INCHIKEY': 'inchikey',
    'KEGG': ['kegg.compound', 'kegg.drug'],
    'KEGG_ID': 'kegg.compound',
    'Labcorp LOINC ID': 'loinc',
    'LM_ID': 'lm',
    'PUBCHEM': 'pubchem.compound',
    'PubChem_CID': 'pubchem.compound',
    'Quest LOINC ID': 'loinc',
    'RefMet_ID': 'rm',
    'SMILES': 'smiles',
    'uniprot': 'uniprotkb'
}

IDNORM = IdentifierNorm(biolink_version="4.2.5")

MAPPING_INPUT_DIR = PROJECT_ROOT_PATH / 'input_data' / 'mapping'
INTERMEDIATE_RESULTS_DIR = PROJECT_ROOT_PATH / 'src' / 'mapping' / 'results_intermediate'


def get_curies(row: pd.Series, id_cols: List[str], array_delimiters: List[str]) -> Tuple[List[str], Dict[str, List[str]]]:
    # First load IDs of different types, handling multiple IDs in one cell as necessary
    local_ids_map = {id_col: [local_id for local_id in re.split(f"[{''.join(array_delimiters)}]", row[id_col])]
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


def annotate_with_curies(dataset_df: pd.DataFrame, id_cols: List[str], dataset_name: str, entity_type: str, version: str, array_delimiters: Optional[List[str]] = None):
    print(f"On dataset {dataset_name}, {entity_type}, {version}")
    array_delimiters = array_delimiters if array_delimiters else [',', ';']  # Default to comma and semicolon

    dataset_df[['curies', 'invalid_ids']] = dataset_df.apply(lambda row: get_curies(row, id_cols, array_delimiters),
                                                             axis=1,
                                                             result_type='expand')
    output_path = INTERMEDIATE_RESULTS_DIR / f"{dataset_name}_{entity_type}_{version}_with_curies.tsv"
    dataset_df.to_csv(output_path, sep='\t', index=False)
    print(dataset_df)


def load_arivale_metabolites_original() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['CAS', 'KEGG', 'HMDB', 'PUBCHEM']
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'metabolomics_metadata.tsv'
    df = pd.read_table(input_path, dtype={'PUBCHEM': str}, skiprows=13)
    return df, id_cols, delimiters

def load_arivale_metabolites_gdprefmet() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'Arivale_GDP_metadata_REFMETANNOT.tsv'
    df = pd.read_table(input_path)
    df[id_cols] = df[id_cols].replace('-', np.nan)
    return df, id_cols, delimiters

def load_arivale_metabolites_lancerefmet() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['CAS', 'KEGG', 'HMDB', 'PUBCHEM', 'INCHIKEY', 'SMILES', 'PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [',', ';']
    # This one is already the source metadata merged with Lance's refmet annotations
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'metabolites_ChemicalAnnotation-Table1.tsv'
    df = pd.read_table(input_path)
    df[id_cols] = df[id_cols].replace('-', np.nan)
    return df, id_cols, delimiters

def load_arivale_proteins_original() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['uniprot', 'gene_id']
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'proteomics_metadata.tsv'
    df = pd.read_table(input_path, skiprows=13)
    return df, id_cols, delimiters

def load_arivale_clinicallabs_original() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['Labcorp LOINC ID', 'Quest LOINC ID']
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'chemistries_metadata.tsv'
    df = pd.read_table(input_path, skiprows=13)
    return df, id_cols, delimiters

def load_arivale_lipids_initial() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['HMDB', 'KEGG']
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'lipids_ChemicalAnnotation-Table1.tsv'
    df = pd.read_table(input_path)
    df[id_cols] = df[id_cols].replace('-', np.nan)
    return df, id_cols, delimiters

def load_arivale_lipids_refmet() -> Tuple[pd.DataFrame, List[str], List[str]]:
    id_cols = ['HMDB', 'KEGG', 'PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [',', ';']

    # Join the source metadata with the refmet annotations file from Lance
    input_path_source = MAPPING_INPUT_DIR / 'arivale' / 'lipids_ChemicalAnnotation-Table1.tsv'
    input_path_refmet = MAPPING_INPUT_DIR / 'arivale' / 'Arivale_lipidomics_metadata_REFMETANNOT.tsv'
    source_df = pd.read_table(input_path_source)
    refmet_df = pd.read_table(input_path_refmet)
    df = pd.merge(source_df, refmet_df, left_on='CHEMICAL_NAME', right_on='Input.name')

    df[id_cols] = df[id_cols].replace('-', np.nan)
    return df, id_cols, delimiters


def main():
    files_to_map = {
        'arivale': {
            'metabolites': {
                'original': load_arivale_metabolites_original,
                'gdprefmet': load_arivale_metabolites_gdprefmet,
                'lancerefmet': load_arivale_metabolites_lancerefmet
            },
            'proteins': {
                'original': load_arivale_proteins_original
            },
            'clinicallabs': {
                'original': load_arivale_clinicallabs_original
            },
            'lipids': {
                'initial': load_arivale_lipids_initial,
                'refmet': load_arivale_lipids_refmet
            }
        }
    }

    for dataset, info in files_to_map.items():
        for entity_type, versions in info.items():
            for version, loader in versions.items():
                df, id_cols, array_delimiters = loader()
                annotate_with_curies(df, id_cols, dataset, entity_type, version, array_delimiters)




if __name__ == "__main__":
    main()
