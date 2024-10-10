import sys
import os.path
import cProfile
import re
from datetime import date
import json
from os import path
from copy import deepcopy
import shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))


from seqteleporter.utils.load_input_params import (validate_input_params, transform_input_excel_sheet_to_text_input,
                                                   show_input_seq_info)
from seqteleporter.post_partition_processor.post_partition_processor import post_partition_processing
from seqteleporter.partitioner.compute_best_partitions import (compute_best_partitions)
from seqteleporter.utils.load_input_params import load_input_params
from seqteleporter.utils.utils import (prepare_0idx_mutations)

print(f'Here is: {SCRIPT_DIR}')
max_partition_number_checked=float('inf')
search_method = 'BFS'


def prepare_compute_best_partitions_params(input_file_path):
    input_params = load_input_params(input_file_path=input_file_path, supress_output=True)
    print(input_params['mutations_1idx'])
    all_mutations_0idx, linked_mutations_0idx = prepare_0idx_mutations(
        input_params['mutations_1idx'], input_params['linked_mutations_1idx']
    )
    params = deepcopy(input_params)
    keys = ['gene_name', 'mutations_1idx', 'linked_mutations_1idx', 'five_prime_dna', 'three_prime_dna',
            'fix_wt_dna_sequence', 'module_plate_format']
    for k in keys:
        params.pop(k, None)
    output_dir = re.sub('[.]txt', f'_{search_method}', input_file_path)
    output_dir = re.sub('__', '_', output_dir)
    output_dir = f'{output_dir}_{date.today()}output'
    params.update(dict(
        mutations_0idx=all_mutations_0idx,
        linked_mutations_0idx=linked_mutations_0idx,
        output_dir=output_dir,
        supress_output=False,
        search_method=search_method,
        sort_by_cost=True,
        partition_search_mode='exhaustive',
        select_top_n_partitions=3,
        max_partition_number_checked=max_partition_number_checked
    ))
    return params


def generate_and_optimize_ready_to_click_modules(text_input_path):
    print('\033[1m=============================================================================================\033[0m')
    print('                               \033[1m RUN STARTED! \033[0m ')
    print('\033[1m=============================================================================================\033[0m')
    validate_input_params(text_input_path)
    show_input_seq_info(text_input_path)
    params = prepare_compute_best_partitions_params(text_input_path)
    print(params['search_method'])
    best_partitions_by_cut_number, best_partitions_by_cut_number_file_path = compute_best_partitions(**params)

    print(f'\n\033[1m'
          f'\n=================================================================================================='
          f'\n                        GENERATE & OPTIMIZE READY-TO-CLICK MODULES '
          f'\n                                     RUN COMPLETED!'
          f'\n=================================================================================================='
          f'\033[0m')


def list_input_files(input_dirs):
    sel_input_files = []
    for input_dir in input_dirs:
        for item in os.listdir(os.path.join(SCRIPT_DIR,  input_dir)):
            full_path = os.path.join(SCRIPT_DIR,  input_dir, item)
            if os.path.isfile(full_path) and re.search('[.]txt$',item):
                sel_input_files.append(full_path)
    return sel_input_files


if __name__=="__main__":
    input_dirs = ['rep0']

    input_files = list_input_files(input_dirs=input_dirs)[3:4]
    for input_file in input_files:
        out_f = f"profile_stats_partitions_{search_method}_{re.sub('[.]txt', '', os.path.basename(input_file))}"
        cProfile.run('generate_and_optimize_ready_to_click_modules(input_file)', path.join(path.dirname(input_file), out_f))
