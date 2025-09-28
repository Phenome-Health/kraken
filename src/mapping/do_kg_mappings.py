import ast
import logging
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Optional, List

import pandas as pd


PROJECT_ROOT_PATH = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT_PATH))

from src.utils.logging_config import setup_logging
from src.utils.kg_io import load_equivalency_mappings

setup_logging()

MAPPING_DIR = PROJECT_ROOT_PATH / 'src' / 'mapping'
INTERMEDIATE_RESULTS_DIR = MAPPING_DIR / 'results_intermediate'
FINAL_RESULTS_DIR = MAPPING_DIR / 'results'


def get_canonical_ids(curies_list, equivalency_map) -> List[str]:
    return list({equivalency_map[curie] for curie in curies_list if curie in equivalency_map})


def get_majority_canonical_id(curies_list, equivalency_map) -> Optional[str]:
    votes = [equivalency_map[curie] for curie in curies_list if curie in equivalency_map]
    if votes:
        most_frequent = Counter(votes).most_common(1)[0][0]
        return most_frequent
    else:
        return None


def main():
    files_to_evaluate = sorted(os.listdir(INTERMEDIATE_RESULTS_DIR))
    print(f"Found {len(files_to_evaluate)} intermediate results files to evaluate.")
    graphs_to_map_to = {
        'kraken': PROJECT_ROOT_PATH / 'artifacts' / 'integrated' / 'kraken_nodes_1.0.1.jsonl',
        'kg2': PROJECT_ROOT_PATH / 'artifacts' / 'harmonized' / 'kg2' / 'nodes.jsonl',
        'spoke': PROJECT_ROOT_PATH / 'artifacts' / 'harmonized' / 'spoke' / 'nodes.jsonl'
    }
    print(f"Will map them to {len(graphs_to_map_to)} graphs: {list(graphs_to_map_to.keys())}")

    for graph_name, nodes_path in graphs_to_map_to.items():
        print(f"On graph {graph_name}")
        equivalency_map = load_equivalency_mappings(nodes_path)
        kg_results_dir = FINAL_RESULTS_DIR / graph_name
        os.makedirs(kg_results_dir, exist_ok=True)

        stats = []
        for file_name in files_to_evaluate:
            print(f"On file {file_name}")
            file_path = INTERMEDIATE_RESULTS_DIR / file_name
            dataset, entity_type, version = file_name.split('_')[:3]
            df = pd.read_table(file_path)
            df.curies = df.curies.apply(ast.literal_eval)
            df.curies_provided = df.curies_provided.apply(ast.literal_eval)
            df.curies_assigned = df.curies_assigned.apply(ast.literal_eval)
            df.invalid_ids = df.invalid_ids.apply(ast.literal_eval)
            df.invalid_ids_provided = df.invalid_ids_provided.apply(ast.literal_eval)
            df.invalid_ids_assigned = df.invalid_ids_assigned.apply(ast.literal_eval)
            df['kg_canonical_ids'] = df.curies.apply(lambda x: get_canonical_ids(x, equivalency_map))
            df['kg_majority_canonical_id'] = df.curies.apply(lambda x: get_majority_canonical_id(x, equivalency_map))
            df['kg_canonical_ids_provided'] = df.curies_provided.apply(lambda x: get_canonical_ids(x, equivalency_map))
            df['kg_canonical_ids_assigned'] = df.curies_assigned.apply(lambda x: get_canonical_ids(x, equivalency_map))

            total_items = len(df)
            has_valid_ids = df.curies.apply(lambda x: len(x) > 0).sum()
            has_valid_ids_provided = df.curies_provided.apply(lambda x: len(x) > 0).sum()
            has_valid_ids_assigned = df.curies_assigned.apply(lambda x: len(x) > 0).sum()
            has_only_provided_ids = has_valid_ids - has_valid_ids_assigned
            has_only_assigned_ids = has_valid_ids - has_valid_ids_provided
            has_both_provided_and_assigned_ids = has_valid_ids - has_only_provided_ids - has_only_assigned_ids
            has_no_ids = ((df.curies.apply(len) == 0) & (df.invalid_ids.apply(len) == 0)).sum()
            has_invalid_ids = df.invalid_ids.apply(lambda x: len(x) > 0).sum()
            has_invalid_ids_provided = df.invalid_ids_provided.apply(lambda x: len(x) > 0).sum()
            has_invalid_ids_assigned = df.invalid_ids_assigned.apply(lambda x: len(x) > 0).sum()
            mapped_to_kg = df.kg_canonical_ids.apply(lambda x: len(x) > 0).sum()
            mapped_to_kg_provided = df.kg_canonical_ids_provided.apply(lambda x: len(x) > 0).sum()
            mapped_to_kg_assigned = df.kg_canonical_ids_assigned.apply(lambda x: len(x) > 0).sum()
            assigned_mappings_verified = df.apply(lambda r: len(set(r.kg_canonical_ids_provided) & set(r.kg_canonical_ids_assigned)) > 0, axis=1).sum()
            has_invalid_ids_and_not_in_kg = ((df.invalid_ids.apply(len) > 0) & (df.kg_canonical_ids.apply(len) == 0)).sum()
            one_to_one_mappings = df.kg_canonical_ids.apply(lambda x: len(x) == 1).sum()
            one_to_many_mappings = df.kg_canonical_ids.apply(lambda x: len(x) > 1).sum()

            file_shortname = f"{dataset}_{entity_type}_{version}"

            headers = ['dataset', 'total_items', 'has_valid_ids', 'has_valid_ids_provided', 'has_valid_ids_assigned',
                       'has_only_provided_ids', 'has_only_assigned_ids', 'has_both_provided_and_assigned_ids',
                       'mapped_to_kg', 'one_to_one_mappings', 'one_to_many_mappings',
                       'mapped_to_kg_provided', 'mapped_to_kg_assigned', 'assigned_mappings_verified',
                       'has_invalid_ids', 'has_invalid_ids_provided', 'has_invalid_ids_assigned',
                       'has_no_ids', 'has_invalid_ids_and_not_in_kg']
            row = [file_shortname, total_items, has_valid_ids, has_valid_ids_provided, has_valid_ids_assigned,
                   has_only_provided_ids, has_only_assigned_ids, has_both_provided_and_assigned_ids,
                   mapped_to_kg, one_to_one_mappings, one_to_many_mappings,
                   mapped_to_kg_provided, mapped_to_kg_assigned, assigned_mappings_verified,
                   has_invalid_ids, has_invalid_ids_provided, has_invalid_ids_assigned,
                   has_no_ids, has_invalid_ids_and_not_in_kg]
            stats.append(row)

            # Record the full results
            df.to_csv(kg_results_dir / f"{file_shortname}_a_full_results.tsv", sep='\t')

            # Record the items that had valid curies but that weren't in the KG, for easy reference
            kg_misses = df[(df.curies.apply(len) > 0) & (df.kg_canonical_ids.apply(len) == 0)]
            kg_misses.to_csv(kg_results_dir / f"{file_shortname}_b_curie_misses.tsv", sep='\t')

            # Record the items that didn't get mapped to the KG, for easy reference
            unmapped = df[df.kg_canonical_ids.apply(len) == 0]
            unmapped.to_csv(kg_results_dir / f"{file_shortname}_c_unmapped.tsv", sep='\t')

            # Record the items that DID map to the KG, for easy reference
            mapped = df[df.kg_canonical_ids.apply(len) > 0]
            mapped.to_csv(kg_results_dir / f"{file_shortname}_d_mapped.tsv", sep='\t')

            # Record the items with invalid IDs, for easy reference
            invalid_ids = df[df.invalid_ids.apply(lambda x: len(x) > 0)]
            invalid_ids.to_csv(kg_results_dir / f"{file_shortname}_e_invalid_ids.tsv", sep='\t')

            # Record the one-to-many items, for easy reference
            one_to_many_items = df[df.kg_canonical_ids.apply(lambda x: len(x) > 1)]
            one_to_many_items.to_csv(kg_results_dir / f"{file_shortname}_f_one_to_many.tsv", sep='\t')


        stats_df = pd.DataFrame(stats, columns=headers)
        summary_stats_path = kg_results_dir / 'a_summary.tsv'
        stats_df.to_csv(summary_stats_path, sep='\t', index=False)

        print(f"Generating visualizations for {graph_name}")
        os.system(f"uv run python {MAPPING_DIR / 'visualize_summary_heatmap.py'} {graph_name.upper()} --input {summary_stats_path} --output {kg_results_dir}")
        os.system(f"uv run python {MAPPING_DIR / 'visualize_stacked_bars_grid.py'} {graph_name.upper()} --input {summary_stats_path} --output {kg_results_dir}")


    # Generate comparative delta barcharts of percent gained/lost with different KGs
    print(f"Generating delta barcharts")
    kraken_summary_path = MAPPING_DIR / 'results' / 'kraken' / 'a_summary.tsv'
    kg2_summary_path = MAPPING_DIR / 'results' / 'kg2' / 'a_summary.tsv'
    spoke_summary_path = MAPPING_DIR / 'results' / 'spoke' / 'a_summary.tsv'
    barchart_run_command = f"uv run python {MAPPING_DIR / 'visualize_delta_barchart.py'}"
    os.system(f"{barchart_run_command} {kg2_summary_path} {kraken_summary_path} --baseline-name KG2 --comparison-name KRAKEN --output {FINAL_RESULTS_DIR}")
    os.system(f"{barchart_run_command} {spoke_summary_path} {kraken_summary_path} --baseline-name SPOKE --comparison-name KRAKEN --output {FINAL_RESULTS_DIR}")
    os.system(f"{barchart_run_command} {kraken_summary_path} {kg2_summary_path} --baseline-name KRAKEN --comparison-name KG2 --output {FINAL_RESULTS_DIR}")
    os.system(f"{barchart_run_command} {kraken_summary_path} {spoke_summary_path} --baseline-name KRAKEN --comparison-name SPOKE --output {FINAL_RESULTS_DIR}")




if __name__ == "__main__":
    main()
