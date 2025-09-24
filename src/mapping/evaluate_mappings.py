import ast
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT_PATH = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT_PATH))


def main():
    input_dir = PROJECT_ROOT_PATH / 'src' / 'mapping'
    output_dir = PROJECT_ROOT_PATH / 'src' / 'mapping'
    files_to_evaluate = ['arivale_metabolites_with_curies.tsv', 'arivale_proteins_with_curies.tsv',
                         'arivale_clinicallabs_with_curies.tsv', 'ukbb_proteins_with_curies.tsv']

    stats = []
    for file_name in files_to_evaluate:
        file_path = input_dir / file_name
        df = pd.read_table(file_path)
        df.curies = df.curies.apply(ast.literal_eval)
        df.invalid_ids = df.invalid_ids.apply(ast.literal_eval)

        file_shortname = ' '.join(file_name.split('_')[:2])

        total_items = len(df)
        has_valid_curies = df.curies.apply(lambda x: len(x) > 0).sum()
        has_no_ids = ((df.curies.apply(len) == 0) & (df.invalid_ids.apply(len) == 0)).sum()
        has_invalid_ids = df.invalid_ids.apply(lambda x: len(x) > 0).sum()

        row = [file_shortname, total_items, has_valid_curies, has_no_ids, has_invalid_ids]
        print(row)



if __name__ == "__main__":
    main()
