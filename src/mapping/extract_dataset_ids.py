import os
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
    'derived_uniprot': 'uniprotkb',
    'gene_id': 'ensembl',
    'HMDB': 'hmdb',
    'HMDB_ID': 'hmdb',
    'INCHI_KEY': 'inchikey',
    'INCHIKEY': 'inchikey',
    'KEGG': ['kegg.compound', 'kegg.drug'],
    'KEGG_ID': 'kegg.compound',
    'Labcorp LOINC ID': 'loinc',
    'LM_ID': 'lm',
    'loinc_code': 'loinc',
    'original_chebi_id': 'chebi',
    'PUBCHEM': 'pubchem.compound',
    'PubChem_CID': 'pubchem.compound',
    'Quest LOINC ID': 'loinc',
    'RefMet_ID': 'rm',
    'SMILES': 'smiles',
    'source_chebi_id': 'chebi',
    'source_hmdb_id': 'hmdb',
    'source_pubchem_id': 'pubchem.compound',
    'uniprot': 'uniprotkb',
    'UniProt': 'uniprotkb'
}

IDNORM = IdentifierNorm(biolink_version="4.2.5")

MAPPING_INPUT_DIR = PROJECT_ROOT_PATH / 'input_data' / 'mapping'
INTERMEDIATE_RESULTS_DIR = PROJECT_ROOT_PATH / 'src' / 'mapping' / 'results_intermediate'
os.makedirs(INTERMEDIATE_RESULTS_DIR, exist_ok=True)


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


def annotate_with_curies(dataset_df: pd.DataFrame,
                         provided_id_cols: List[str],
                         assigned_id_cols: List[str],
                         dataset_name: str,
                         entity_type: str,
                         version: str,
                         array_delimiters: List[str]):
    print(f"On dataset {dataset_name}, {entity_type}, {version}")
    id_cols = provided_id_cols + assigned_id_cols

    # First get curies for ALL id columns
    dataset_df[['curies', 'invalid_ids']] = dataset_df.apply(lambda row: get_curies(row, id_cols, array_delimiters),
                                                             axis=1,
                                                             result_type='expand')

    # Then get them for two subsets of the id columns (those provided in the dataset and those we assigned somehow)
    dataset_df[['curies_provided', 'invalid_ids_provided']] = dataset_df.apply(lambda row: get_curies(row, provided_id_cols, array_delimiters),
                                                                               axis=1,
                                                                               result_type='expand')

    # Then get them for two subsets of the id columns (those provided in the dataset and those we assigned somehow)
    dataset_df[['curies_assigned', 'invalid_ids_assigned']] = dataset_df.apply(lambda row: get_curies(row, assigned_id_cols, array_delimiters),
                                                                               axis=1,
                                                                               result_type='expand')

    # Save the intermediate results file
    output_path = INTERMEDIATE_RESULTS_DIR / f"{dataset_name}_{entity_type}_{version}_with_curies.tsv"
    dataset_df.to_csv(output_path, sep='\t', index=False)
    print(dataset_df)


def load_arivale_metabolites_original() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['CAS', 'KEGG', 'HMDB', 'PUBCHEM']
    assigned_id_cols = []
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'metabolomics_metadata.tsv'
    df = pd.read_table(input_path, dtype={'PUBCHEM': str}, skiprows=13)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_arivale_metabolites_gdprefmet() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'Arivale_GDP_metadata_REFMETANNOT.tsv'
    df = pd.read_table(input_path)
    df[assigned_id_cols] = df[assigned_id_cols].replace('-', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_arivale_metabolites_lancerefmet() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['CAS', 'KEGG', 'HMDB', 'PUBCHEM', 'INCHIKEY', 'SMILES']
    assigned_id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [',', ';']
    # This one is already the source metadata merged with Lance's refmet annotations
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'metabolites_ChemicalAnnotation-Table1.tsv'
    df = pd.read_table(input_path)
    df[assigned_id_cols] = df[assigned_id_cols].replace('-', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_arivale_proteins_original() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['uniprot', 'gene_id']
    assigned_id_cols = []
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'proteomics_metadata.tsv'
    df = pd.read_table(input_path, skiprows=13)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_arivale_clinicallabs_original() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['Labcorp LOINC ID', 'Quest LOINC ID']
    assigned_id_cols = []
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'chemistries_metadata.tsv'
    df = pd.read_table(input_path, skiprows=13)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_arivale_lipids_initial() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['HMDB', 'KEGG']
    assigned_id_cols = []
    delimiters = [',', ';']
    input_path = MAPPING_INPUT_DIR / 'arivale' / 'lipids_ChemicalAnnotation-Table1.tsv'
    df = pd.read_table(input_path)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_arivale_lipids_refmet() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['HMDB', 'KEGG']
    assigned_id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [',', ';']

    # Join the source metadata with the refmet annotations file from Lance
    input_path_source = MAPPING_INPUT_DIR / 'arivale' / 'lipids_ChemicalAnnotation-Table1.tsv'
    input_path_refmet = MAPPING_INPUT_DIR / 'arivale' / 'Arivale_lipidomics_metadata_REFMETANNOT.tsv'
    source_df = pd.read_table(input_path_source)
    refmet_df = pd.read_table(input_path_refmet)
    df = pd.merge(source_df, refmet_df, left_on='CHEMICAL_NAME', right_on='Input.name')

    df[assigned_id_cols] = df[assigned_id_cols].replace('-', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_ukbb_proteins_original() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = ['UniProt']
    assigned_id_cols = []
    delimiters = ['_']
    input_path = MAPPING_INPUT_DIR / 'ukbb' / 'UKBB_Protein_Meta.tsv'
    df = pd.read_table(input_path)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_ukbb_clinicallabs_filtered() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['loinc_code']
    delimiters = ['_']

    # The biomapper results for this one contain qc fields as well; filter those out
    biomapper_results_path = MAPPING_INPUT_DIR / 'ukbb' / 'ukbb_chemistry_COMPLETE.tsv'
    fieldnames_path = MAPPING_INPUT_DIR / 'ukbb' / 'clinicallab_fieldnames.tsv'  # Copied from web (per Lance)
    fieldnames_df = pd.read_table(fieldnames_path)
    field_names = set(fieldnames_df['Description'].unique())
    print(f"Field names to filter to: {field_names}")
    df = pd.read_table(biomapper_results_path)
    df = df[df['field_name'].isin(field_names)]

    df[assigned_id_cols] = df[assigned_id_cols].replace('NO_MATCH', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_ukbb_metabolites_biomapper() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['source_chebi_id', 'source_hmdb_id', 'source_pubchem_id']  # TODO: verify these are assigned?
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'ukbb' / 'ukbb_metabolites_COMPLETE.tsv'
    df = pd.read_table(input_path, dtype={'source_pubchem_id': str})
    df.source_pubchem_id = df.source_pubchem_id.str.removesuffix('.0')
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_israeli10k_lipids_refmet() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'israeli10k' / 'israeli10k_lipidomics_metadata_REFMETANNOT.tsv'
    df = pd.read_table(input_path)
    df[assigned_id_cols] = df[assigned_id_cols].replace('-', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_israeli10k_lipids_website() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['PubChem_CID', 'ChEBI_ID', 'HMDB_ID', 'LM_ID', 'KEGG_ID', 'INCHI_KEY', 'RefMet_ID']
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'israeli10k' / 'Israeli10k_website_lipidomics_metadata_REFMETANNOT.tsv'
    df = pd.read_table(input_path)
    df[assigned_id_cols] = df[assigned_id_cols].replace('-', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_israeli10k_metabolites_biomapper() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['original_chebi_id']  # TODO: Verify this is assigned?
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'israeli10k' / 'israeli10k_metabolites_COMPLETE.tsv'
    df = pd.read_table(input_path)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_israeli10k_proteins_biomapper() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['derived_uniprot']
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'israeli10k' / 'israeli10k_nightingale_proteins_mapped.tsv'
    # NOTE: I manually collapsed the one derived measure into one in this tsv; has both uniprot ids
    df = pd.read_table(input_path)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_israeli10k_clinicallabs_biomapper() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['loinc_code']
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'israeli10k' / 'israeli10k_chemistry_loinc_COMPLETE.tsv'
    df = pd.read_table(input_path)
    return df, provided_id_cols, assigned_id_cols, delimiters

def load_israeli10k_clinicallabs_biomapper2() -> Tuple[pd.DataFrame, List[str], List[str], List[str]]:
    provided_id_cols = []
    assigned_id_cols = ['loinc_code']
    delimiters = [';']
    input_path = MAPPING_INPUT_DIR / 'israeli10k' / 'israeli10k_chemistry_loinc_COMPLETE_2.tsv'
    df = pd.read_table(input_path)
    df[assigned_id_cols] = df[assigned_id_cols].replace('NO_MATCH', np.nan)
    return df, provided_id_cols, assigned_id_cols, delimiters


def main():
    files_to_map = {
        'arivale': {
            'metabolites': {
                'v1_original': load_arivale_metabolites_original,
                'v2_gdprefmet': load_arivale_metabolites_gdprefmet,
                'v3_lancerefmet': load_arivale_metabolites_lancerefmet
            },
            'proteins': {
                'v1_original': load_arivale_proteins_original
            },
            'clinicallabs': {
                'v1_original': load_arivale_clinicallabs_original
            },
            'lipids': {
                'v1_initial': load_arivale_lipids_initial,
                'v2_refmet': load_arivale_lipids_refmet
            }
        },
        'ukbb': {
            'proteins': {
                'v1_original': load_ukbb_proteins_original
            },
            'clinicallabs': {
                'v1_filtered': load_ukbb_clinicallabs_filtered
            },
            'metabolites': {
                'v1_biomapper': load_ukbb_metabolites_biomapper
            }
        },
        'israeli10k': {
            'lipids': {
                'v1_refmet': load_israeli10k_lipids_refmet,
                'v2_websiterefmet': load_israeli10k_lipids_website
            },
            'metabolites': {
                'v1_biomapper': load_israeli10k_metabolites_biomapper
            },
            'proteins': {
                'v1_biomapper': load_israeli10k_proteins_biomapper
            },
            'clinicallabs': {
                'v1_biomapper': load_israeli10k_clinicallabs_biomapper,
                'v2_biomapper': load_israeli10k_clinicallabs_biomapper2
            }
        }
    }

    # Loop through and extract ids from all the above-defined files
    for dataset, info in files_to_map.items():
        for entity_type, versions in info.items():
            for version, loader in versions.items():
                df, provided_id_cols, assigned_id_cols, array_delimiters = loader()
                annotate_with_curies(df, provided_id_cols, assigned_id_cols, dataset, entity_type, version, array_delimiters)



if __name__ == "__main__":
    main()
