import ast
import logging
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT_PATH = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT_PATH))

from src.utils.logging_config import setup_logging
from src.utils.kg_io import load_equivalency_mappings

setup_logging()


def get_canonical_ids(curies_list, equivalency_map):
    return list({equivalency_map[curie] for curie in curies_list if curie in equivalency_map})


def main():
    input_dir = PROJECT_ROOT_PATH / 'src' / 'mapping'
    output_dir = PROJECT_ROOT_PATH / 'src' / 'mapping' / 'results'
    files_to_evaluate = ['arivale_metabolites_with_curies.tsv', 'arivale_proteins_with_curies.tsv',
                         'arivale_clinicallabs_with_curies.tsv', 'ukbb_proteins_with_curies.tsv']


    equivalency_map = load_equivalency_mappings(PROJECT_ROOT_PATH / 'artifacts' / 'integrated' / 'kraken_nodes_1.0.1.jsonl')


    headers = ['dataset', 'total_items', 'has_valid_ids', 'in_kraken', 'no_provided_ids', 'has_invalid_ids', 'has_invalid_ids_and_not_in_kraken']
    stats = []
    for file_name in files_to_evaluate:
        file_path = input_dir / file_name
        df = pd.read_table(file_path)
        df.curies = df.curies.apply(ast.literal_eval)
        df.invalid_ids = df.invalid_ids.apply(ast.literal_eval)
        df['kraken_canonical_ids'] = df.curies.apply(lambda x: get_canonical_ids(x, equivalency_map))

        file_shortname = '_'.join(file_name.split('_')[:2])

        total_items = len(df)
        has_valid_curies = df.curies.apply(lambda x: len(x) > 0).sum()
        has_no_ids = ((df.curies.apply(len) == 0) & (df.invalid_ids.apply(len) == 0)).sum()
        has_invalid_ids = df.invalid_ids.apply(lambda x: len(x) > 0).sum()
        in_kraken = df.kraken_canonical_ids.apply(lambda x: len(x) > 0).sum()
        has_invalid_ids_and_not_in_kraken = ((df.invalid_ids.apply(len) > 0) & (df.kraken_canonical_ids.apply(len) == 0)).sum()

        row = [file_shortname, total_items, has_valid_curies, in_kraken, has_no_ids, has_invalid_ids, has_invalid_ids_and_not_in_kraken]
        stats.append(row)

        # Record the full results
        df.to_csv(output_dir / f"{file_shortname}_a_full_results.tsv", sep='\t')

        # Record the items that had valid curies but that weren't in kraken, for easy reference
        kraken_misses = df[(df.curies.apply(len) > 0) & (df.kraken_canonical_ids.apply(len) == 0)]
        kraken_misses.to_csv(output_dir / f"{file_shortname}_b_curie_misses.tsv", sep='\t')

        # Record the items that didn't get mapped to kraken, for easy reference
        unmapped = df[df.kraken_canonical_ids.apply(len) == 0]
        unmapped.to_csv(output_dir / f"{file_shortname}_c_unmapped.tsv", sep='\t')

        # Record the items that DID map to kraken, for easy reference
        mapped = df[df.kraken_canonical_ids.apply(len) > 0]
        mapped.to_csv(output_dir / f"{file_shortname}_d_mapped.tsv", sep='\t')

        # Record the items with invalid IDs, for easy reference
        invalid_ids = df[df.invalid_ids.apply(lambda x: len(x) > 0)]
        invalid_ids.to_csv(output_dir / f"{file_shortname}_e_invalid_ids.tsv", sep='\t')


    stats_df = pd.DataFrame(stats, columns=headers)
    stats_df.to_csv(output_dir / 'summary.tsv', sep='\t', index=False)



if __name__ == "__main__":
    main()
