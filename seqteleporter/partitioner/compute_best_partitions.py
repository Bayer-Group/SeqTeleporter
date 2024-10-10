import json
import re
import time
from datetime import date
from copy import deepcopy
from itertools import chain
from os import makedirs, listdir, remove
from os.path import join, exists
from typing import Tuple, Union, List

import pandas as pd

from seqteleporter.config import ENZYME_INFO, PARTITION_SEARCH_MODES
from seqteleporter.partition_property_finder.fusion_sites_finder import concat_sel_fusion_sites_to_fragments
from seqteleporter.partition_property_finder.partition_property_finder import find_partition_property
from seqteleporter.partitioner.partitioner import find_cuttable_positions, partitioner
from seqteleporter.utils.load_input_params import load_input_params
from seqteleporter.utils.utils import (zero_indexing_to_one_indexing, compute_lib_complexity,
                                       is_valid_fusion_site_set, is_dna, is_aa, pretty_fragments_expression,
                                       annotate_mutations_in_fragments, prepare_0idx_mutations)


def prepare_output_dirs(output_dir):
    """create output dirs if not exist"""

    log_dir = join(output_dir, 'logs')
    result_dir = join(output_dir, 'results')
    dirs_to_create = [output_dir, log_dir, result_dir]
    for directory in dirs_to_create:
        if not exists(directory):
            makedirs(directory)

    # clean log directory before starting the run
    for f in listdir(log_dir):
        if not re.search('ipynb_checkpoints', f):
            remove(join(log_dir, f))
    return log_dir, result_dir


def prepare_compute_best_partitions_params(input_file_path: str) -> dict:
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
    output_dir = re.sub('[.]txt', '', input_file_path)
    output_dir = re.sub('__', '_', output_dir)
    output_dir = f'{output_dir}_{date.today()}output'
    params.update(dict(
        mutations_0idx=all_mutations_0idx,
        linked_mutations_0idx=linked_mutations_0idx,
        output_dir=output_dir,
        supress_output=False,
        search_method="BFS",
        sort_by_cost=True,
        partition_search_mode='dist_mut_1',
        select_top_n_partitions=3,
        max_partition_number_checked=100000
    ))
    return params


"""Define function - Computes the best possible partitions of a given sequence based on a range of criteria including 
the number of cuts, minimum fragment length, maximum cost, maximum unevenness in fragment lengths, and minimum 
ligation fidelity. This function iterates over a specified range of cut numbers and generates partitions that satisfy 
the given constraints, utilizing experimental data for fusion site information."""


def write_compute_best_partitions_log_header(
        s: str,
        mutations_0idx: Union[list, None],
        linked_mutations_0idx: Union[list, None],
        fusion_sites_used_by_backbone: Tuple[str, ...],
        min_aa_length: int,
        max_cost: int,
        max_unevenness: float,
        min_ligation_fidelity: float,
        satisfaction_fidelity: float,
        search_method: str,
        host: str,
        sort_by_cost: bool,
        compute_best_partitions_log_file_path: str
) -> None:

    with open(compute_best_partitions_log_file_path, 'a') as fd:
        fd.write(f'\n Sequence: {s}'
                 f'\n Mutations(0-indexed): {mutations_0idx}'
                 f'\n Linked mutations(0-indexed): {linked_mutations_0idx}'
                 f'\n Fusion sites used by backbone {fusion_sites_used_by_backbone}'
                 f'\n Min length: {min_aa_length}'
                 f'\n Max cost: {max_cost}'
                 f'\n Max unevenness: {max_unevenness}'
                 f'\n Min ligation fidelity: {min_ligation_fidelity}'
                 f'\n Satisfaction fidelity: {satisfaction_fidelity}'
                 f'\n Search method: {search_method}'
                 f'\n Host: {host}'
                 f'\n Sort by cost: {sort_by_cost}')


def write_compute_best_partitions_log_body(
        compute_best_partitions_log_file_path: str,
        number_of_cuts: int,
        elapsed_time_number_of_cuts: float,
        num_of_checked_partitions: int,
        num_of_checked_unique_partitions: int,
        hard_constraint_violations: dict,
        select_top_n_partitions: int,
        sel_partitions: List[dict],
        mutations_0idx: Union[list, None],
        linked_mutations_0idx: Union[list, None],
        supress_output: bool
) -> None:

    with open(compute_best_partitions_log_file_path, 'a') as fd:
        compute_best_partitions_log_header = \
            (f'\n--------------------------------------------------------------------------------------------------'
             f'\n\033[1m Results For {number_of_cuts + 1} Fragments: \033[0m'
             f'\n   Number of cuts: {number_of_cuts} '
             f'\n   Elapsed time: {elapsed_time_number_of_cuts} '
             f'\n   Number of partitions checked: {num_of_checked_partitions}'
             f'\n   Number of unique partitions checked: {num_of_checked_unique_partitions}'
             f'\n   Number of unique partitions violating hard constraints: '
             f'{sum([v for v in hard_constraint_violations.values()])}'
             f'\n   {pd.DataFrame(hard_constraint_violations, index=["Partition Count"]).transpose()}'
             f'\n   The top {select_top_n_partitions} partitions:')

        fd.write(compute_best_partitions_log_header)
        if not supress_output:
            print(compute_best_partitions_log_header)
        for idx, d in enumerate(sel_partitions):
            pretty_fragments = annotate_mutations_in_fragments(fragments=d['fragments'],
                                                               mutations_0idx=mutations_0idx,
                                                               linked_mutations_0idx=linked_mutations_0idx)
            pretty_fragments_txt = ']\n['.join(pretty_fragments)
            compute_best_partitions_log_txt = \
                (
                    f"\n \033[1mPartition {idx + 1}:\033[0m"
                    f"\n \033[1mfragments = \n[{pretty_fragments_txt}]\033[0m"
                    f"\n fusion_sites = {d['fusion_sites']};"
                    f"\n partition(0idx) = {d['partition']}; fidelity={d['ligation_fidelity']}; "
                    f"unevenness = {d['fragment_length_unevenness']}; cost={d['cost']};"
                    f"\n aa + dna fusion sites = {d['expression']}"
                )
            fd.write(compute_best_partitions_log_txt)
            if not supress_output:
                print(compute_best_partitions_log_txt)


def validate_inputs(s: str, fusion_sites_used_by_backbone: Tuple[str, ...]) -> None:
    if not is_aa(s):
        raise ValueError(f"The provided input sequence is not a valid amino acid sequence!")

    if sum([is_dna(fs) for fs in fusion_sites_used_by_backbone]) != len(fusion_sites_used_by_backbone):
        raise ValueError(f"The fusion sites used by backbone are invalid! "
                         f"Some of the fusion sites are not DNA sequence.")

    if not is_valid_fusion_site_set(fusion_sites_used_by_backbone):
        raise ValueError(f"The fusion sites used by backbone are invalid! "
                         f"Some of the fusion sites (or its reverse complement) used by backbone is duplicated."
                         f"Please try again with valid fusion sites.")


def pick_top_n_partitions(res_per_count: List[dict], select_top_n_partitions: int, sort_by_cost: bool) -> List[dict]:
    if sort_by_cost:
        sorted_partitions = sorted(res_per_count, key=lambda x: (
            x['cost'], -x['ligation_fidelity'], x['fragment_length_unevenness'], x['partition']))
    else:
        sorted_partitions = sorted(res_per_count, key=lambda x: (
            -x['ligation_fidelity'], x['cost'], x['fragment_length_unevenness'], x['partition']))

    seen_cut_positions = []
    sel_partitions = []
    for partition in sorted_partitions:
        if partition['partition'] not in seen_cut_positions:
            sel_partitions.append(partition)
            if len(sel_partitions) == select_top_n_partitions:
                break
        seen_cut_positions.append(partition['partition'])

    return sel_partitions


def compute_best_partitions(s: str, mutations_0idx: Union[list, None], linked_mutations_0idx: Union[list, None],
                            cut_number_range: Tuple[int, int], fidelity_data_path: str,
                            fusion_sites_used_by_backbone: Tuple[str, ...], min_aa_length: int,
                            max_cost: int, max_unevenness: float,
                            min_ligation_fidelity: float, satisfaction_fidelity: float, output_dir: str,
                            supress_output: bool, search_method: str, codon_usage_table_path: str, host: str,
                            sort_by_cost: bool, enzyme: str, allowed_cut_positions_1idx: list,
                            partition_search_mode: str, select_top_n_partitions: int,
                            cost_per_nt: float, provider_min_frag_len: int, provider_max_frag_len: int,
                            max_partition_number_checked: int) -> tuple:

    print('\033[1m=============================================================================================\033[0m')
    print('                               \033[1m IDENTIFYING OPTIMAL PARTITIONS \033[0m ')
    print('\033[1m=============================================================================================\033[0m')

    saved_args = {**locals()}
    cost_scan = PARTITION_SEARCH_MODES[partition_search_mode]['cost_scan']
    pre_distribute_mutations = PARTITION_SEARCH_MODES[partition_search_mode]['pre_distribute_mutations']
    one_dist = PARTITION_SEARCH_MODES[partition_search_mode]['one_dist']

    st_total = time.time()
    # first validate the inputs
    validate_inputs(s, fusion_sites_used_by_backbone)

    log_dir, result_dir = prepare_output_dirs(output_dir)
    compute_best_partitions_log_file_path = join(log_dir, 'compute_best_partitions_log.txt')

    write_compute_best_partitions_log_header(s, mutations_0idx, linked_mutations_0idx, fusion_sites_used_by_backbone,
                                             min_aa_length, max_cost, max_unevenness, min_ligation_fidelity,
                                             satisfaction_fidelity, search_method, host, sort_by_cost,
                                             compute_best_partitions_log_file_path)

    cuttable_positions = find_cuttable_positions(s, mutations_0idx, linked_mutations_0idx, min_aa_length,
                                                 provider_max_frag_len, enzyme, allowed_cut_positions_1idx,
                                                 ENZYME_INFO)
    # some constants
    select_top_n_partitions = 3
    res_all = []
    for number_of_cuts in range(cut_number_range[0], cut_number_range[1]):
        st_number_of_cuts = time.time()
        res_per_count: list = []
        max_cost_scan = len(s) * 3 * 0.06
        num_of_checked_partitions = 0
        num_of_checked_unique_partitions = 0
        hard_constraint_violations: dict = {}
        loop_flag = True
        partitions_list = partitioner(s, cuttable_positions, number_of_cuts, mutations_0idx, linked_mutations_0idx,
                                      pre_distribute_mutations, one_dist)
        if isinstance(partitions_list, str):
            hard_constraint_violations.update({partitions_list: 1})
        else:
            while len(res_per_count) == 0 and max_cost_scan <= max_cost and loop_flag:
                # make partitions. partitions have to be generated for every loop, because the output of partitioner is a
                # list of generator-like itertools, which can only be used once.
                partitions_list = partitioner(
                    s,
                    cuttable_positions,
                    number_of_cuts,
                    mutations_0idx,
                    linked_mutations_0idx,
                    pre_distribute_mutations,
                    one_dist
                )

                params = dict(
                    s=s,
                    mutations_0idx=mutations_0idx,
                    linked_mutations_0idx=linked_mutations_0idx,
                    partitions_list=partitions_list,
                    fidelity_data_path=fidelity_data_path,
                    fusion_sites_used_by_backbone=fusion_sites_used_by_backbone,
                    min_aa_length=min_aa_length,
                    max_cost=max_cost,
                    max_unevenness=max_unevenness,
                    min_ligation_fidelity=min_ligation_fidelity,
                    satisfaction_fidelity=satisfaction_fidelity,
                    search_method=search_method,
                    enzyme=enzyme,
                    analyze_cause_of_no_valid_partition=True,
                    cost_per_nt=cost_per_nt,
                    provider_min_frag_len=provider_min_frag_len,
                    provider_max_frag_len=provider_max_frag_len
                )
                if cost_scan:
                    params.update({'max_cost': max_cost_scan})
                else:
                    loop_flag = False

                partition_properties_generator = find_partition_property(**params)
                max_cost_scan += max_cost // 20
                begin, end = [0], [len(s)]
                num_of_checked_unique_partitions = 0
                hard_constraint_violations: dict = {}
                for res in partition_properties_generator:
                    if num_of_checked_unique_partitions >= max_partition_number_checked:
                        break
                    num_of_checked_unique_partitions += 1
                    num_of_checked_partitions += 1
                    if isinstance(res, str):
                        if res not in hard_constraint_violations:
                            hard_constraint_violations.update({res: 1})
                        else:
                            hard_constraint_violations[res] += 1
                    else:
                        lst = [s[sl] for sl in map(slice, chain(begin, res['partition']), chain(res['partition'], end))]
                        res.update({'fragments': lst})
                        res_per_count.append(res)
        if len(res_per_count) == 0:
            sel_partitions = []
        else:
            sel_partitions = pick_top_n_partitions(res_per_count, select_top_n_partitions, sort_by_cost)

            # concat fusion sites back to amino acid after selecting best partitions to avoid running this process in
            # each and every loop of "partitioner()"
            for partition in sel_partitions:
                fragment_with_fusion_sites = concat_sel_fusion_sites_to_fragments(
                    fragments=partition["fragments"],
                    fusion_sites=partition["fusion_sites"],
                    sel_junction_dna_map_fusion_sites=partition['sel_junction_dna_map_fusion_sites'],
                    codon_usage_table_path=codon_usage_table_path
                )
                # create easy-to-read expression of fragments and fusion sites
                if len(partition['partition']) == 0:
                    expr = s
                else:
                    expr = pretty_fragments_expression(
                        fragments=partition["fragments"],
                        fragment_with_fusion_sites=fragment_with_fusion_sites,
                        fusion_site_len=int(ENZYME_INFO[enzyme]['fusion_site_length'])
                    )
                partition.update({
                    'expression': expr,
                    'fragment_with_fusion_sites': fragment_with_fusion_sites
                })
                partition.pop('sel_junction_dna_map_fusion_sites')
                partition.pop("hard_constraint_violation")
        et_number_of_cuts = time.time()
        elapsed_time_number_of_cuts = round(et_number_of_cuts - st_number_of_cuts, 2)
        res_per_cut = {
            'number_of_cuts': number_of_cuts, 'elapsed_time': elapsed_time_number_of_cuts,
            'number_of_partitions_checked': num_of_checked_partitions, 'sel_partitions': sel_partitions,
            'num_of_checked_unique_partitions': num_of_checked_unique_partitions,
            'hard_constraint_violations': hard_constraint_violations
        }
        res_all.append(res_per_cut)

        result_per_cut_path = join(log_dir, f'{partition_search_mode}_{number_of_cuts+1}frags.json')
        with open(result_per_cut_path, 'w') as fp:
            json.dump(res_per_cut, fp)

        write_compute_best_partitions_log_body(compute_best_partitions_log_file_path, number_of_cuts,
                                               elapsed_time_number_of_cuts, num_of_checked_partitions,
                                               num_of_checked_unique_partitions, hard_constraint_violations,
                                               select_top_n_partitions, sel_partitions, mutations_0idx,
                                               linked_mutations_0idx, supress_output)

    diversity = compute_lib_complexity(mutations=mutations_0idx, linked_mutations=linked_mutations_0idx)
    et_total = time.time()
    total_elapsed_time = round(et_total - st_total, 2)
    best_partitions_by_cut_number = {'sequence': s,
                                     'mutations': zero_indexing_to_one_indexing(mutations_0idx),
                                     'total_elapsed_time': total_elapsed_time,
                                     'diversity': diversity,
                                     'best_partitions_by_cut_number': res_all}

    output_path_to_return = join(log_dir, f'best_partitions_by_cut_number_{partition_search_mode}.json')
    for item, name_ in zip([best_partitions_by_cut_number, saved_args], ['best_partitions_by_cut_number', 'params']):
        output_path = join(log_dir, f'{name_}_{partition_search_mode}.json')
        with open(output_path, 'w') as fp:
            json.dump(item, fp)
            print('--------------------------------------------------------------------------------------------------')
            print(f'\033[1m Result "{name_}" is exported at:\n {output_path}\033[0m')

    print(f'\n\033[1m'
          f'\n==========================================================================================='
          f'\n                       SUCCESSFULLY IDENTIFIED OPTIMAL PARTITIONS!                           '
          f'\n                       Total elapsed time: {best_partitions_by_cut_number["total_elapsed_time"]} seconds'
          f'\n===========================================================================================\033[0m')
    return best_partitions_by_cut_number, output_path_to_return


"""Define function - Identifies and returns the lowest cost partition from a collection of partitioning results 
generated by the compute_best_partitions function. This function iterates through the best partitioning result for 
each specified number of cuts, and selects the partition with the lowest cost."""


def get_lowest_cost_from_best_partitions(best_partitions_by_cut_number, supress_output=False):

    res_all = best_partitions_by_cut_number["best_partitions_by_cut_number"]
    sel_cost = float('inf')
    sel_partition = {}
    not_found = True
    for res_per_cut_number in res_all:
        sel_partitions = res_per_cut_number['sel_partitions']
        if len(sel_partitions) > 0:
            not_found = False
            cost = res_per_cut_number['sel_partitions'][0]['cost']
            if cost < sel_cost:
                sel_cost = cost
                sel_partition = res_per_cut_number.copy()
                sel_partition.pop('sel_partitions')
                sel_partition.update(res_per_cut_number['sel_partitions'][0])

    if not not_found:
        sel_partition.update({'sequence': best_partitions_by_cut_number["sequence"],
                              'mutations': best_partitions_by_cut_number["mutations"],
                              'diversity': best_partitions_by_cut_number["diversity"]})

    if not supress_output and not not_found:
        print('--------------------------------------------------------------------------------------------------')
        print('\033[1m The partition with the lowest cost out of all fragment numbers: \033[0m')
        print('--------------------------------------------------------------------------------------------------')
        omit = ['hard_constraint_violations_analysis_when_no_results_are_found',
                'fragment_with_fusion_sites', 'sequence', 'mutations']

        for k, v in sel_partition.items():
            if k not in omit:
                print(k.capitalize(), ': ', v)

    return sel_partition
