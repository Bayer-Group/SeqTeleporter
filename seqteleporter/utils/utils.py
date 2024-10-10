"""Utils Functions"""
from os import path, listdir
import os
import itertools
import pandas as pd
import numpy as np
import math
import re
from typing import Union, Dict, Tuple, Optional, List, Any, Generator
from Bio.Restriction import *
import python_codon_tables as pct
import heapq
from functools import cache

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'


def find_point_mutations(wt_seq: str, mut_seq: str) -> list[str]:
    """Find point mutations"""
    if len(wt_seq) != len(mut_seq):
        raise ValueError(f'The length of the wild type seq is different from the length of the mutant seq!\n'
                         f'Wild type seq length: {len(wt_seq)},\n'
                         f'Mutant seq length: {len(mut_seq)},\n'
                         f'Wild type seq:\n'
                         f'{wt_seq}\n'
                         f'Wild type seq:\n'
                         f'{mut_seq}\n')
    muts = []
    for i in range(0, len(wt_seq)):
        if wt_seq[i] != mut_seq[i]:
            muts.append("".join([wt_seq[i], str(i + 1), mut_seq[i]]))
    return muts


def group_mutations_when_too_many_mutations(mutations: list, num_of_group: int) -> Dict[int, list]:
    new_muts = {}
    for idx, mut in enumerate(mutations):
        if idx % num_of_group not in new_muts:
            new_muts.update({idx % num_of_group: [mut]})
        else:
            new_muts[idx % num_of_group].append(mut)
    return new_muts


def annotate_sequence_mutations(s: str, mutations_0idx: list, linked_mutations_0idx: list) -> dict:
    # https://stackoverflow.com/questions/4842424/list-of-ansi-color-escape-sequences
    ansi_colors = list(range(92, 97)) + list(range(32, 37))
    annotated = {}
    linked_mutation_positions = [mut[1] for muts in linked_mutations_0idx for mut in
                                 muts] if linked_mutations_0idx else []
    for idx, aa in enumerate(s):
        if mutations_0idx:
            if idx in [mut['position'] for mut in mutations_0idx] and idx not in linked_mutation_positions:
                aa = "".join(["\033[31m", aa, "\033[39m"])
        if linked_mutations_0idx:
            for j, muts in enumerate(linked_mutations_0idx):
                if idx in [mut[1] for mut in muts]:
                    annotate_txt = "".join(["\033[", str(ansi_colors[j]), "m"])
                    aa = "".join([annotate_txt, aa, "\033[39m"])
        annotated.update({idx: aa})
    return annotated


def annotate_mutations_in_fragments(fragments: list, mutations_0idx: List[Any],
                                    linked_mutations_0idx: List[Any]) -> List:

    start_positions = [sum([len(f) for f in fragments[0:i]]) for i in range(0, len(fragments) + 1)]
    s = "".join(fragments)
    annotated_seq = annotate_sequence_mutations(s, mutations_0idx, linked_mutations_0idx)
    annotated_frags = []
    for idx, start in enumerate(start_positions):
        if idx < len(start_positions) - 1:
            annotated_f = "".join([annotated_seq[i] for i in range(start, start_positions[idx + 1])])
            annotated_frags.append(annotated_f)
    return annotated_frags


def pretty_fragments_expression(fragments: list, fragment_with_fusion_sites: dict, fusion_site_len: int, ) -> str:
    cumulative_ident_count = 0
    frag_exprs = []
    for f in fragments:
        frag_expr = "".join([
            fragment_with_fusion_sites[f].get("n_term_dna", ""),
            '[',
            fragment_with_fusion_sites[f].get("middle_aa"),
            ']',
            fragment_with_fusion_sites[f].get("c_term_dna", "")
        ])
        frag_exprs.append("".join(["\n", " " * cumulative_ident_count, frag_expr]))
        ident_count = len(frag_expr) - fusion_site_len - 3  # 3 for these symbols: <|>
        cumulative_ident_count = cumulative_ident_count + ident_count
    return "".join(frag_exprs)


def get_available_resources():
    hosts = []
    for f in listdir(path.join(path.dirname(path.dirname(__file__)), 'data', 'codon_usage')):
        if re.match(".*csv$", f):
            hosts.append(re.sub('[.]csv', '', f))

    assembly_conds = []
    for f in listdir(path.join(path.dirname(path.dirname(__file__)), 'data', 'neb_fidelity_data')):
        if re.match(".*xlsx$", f):
            assembly_conds.append(f)
    return {'codon_usage': hosts, 'fidelity_file': assembly_conds}


def print_available_resources():
    print('\033[1m===========================================================================================\033[0m')
    print(f'                                           AVAILABLE HOSTS                                              ')
    print('\033[1m===========================================================================================\033[0m')
    # PRINT THE LIST OF NAMES OF ALL AVAILABLE TABLES
    print(*pct.available_codon_tables_names, sep='\n')

    print('\033[1m===========================================================================================\033[0m')
    print(f'                                  AVAILABLE ASSEMBLY CONDITIONS                                         ')
    print('\033[1m===========================================================================================\033[0m')
    for f in listdir(path.join(path.dirname(path.dirname(__file__)), 'data', 'neb_fidelity_data')):
        if re.match(".*xlsx$", f):
            print(re.sub('FileS[0-9][0-9]_|[.]xlsx', '', f))


def validate_codon_table(codon_table: dict) -> bool:
    # List of single-letter codes for the 20 standard amino acids
    if ''.join(codon_table.keys()) != AMINO_ACIDS:
        raise ValueError('Provided codon table does not contain all 20 amino acids.')
    return True


def validate_enzyme_and_enzyme_info(
        enzyme: str,
        enyzme_info: Dict[str, Union[Dict[str, Union[str, int]], Dict[str, Union[str, int]]]]
) -> bool:
    if enzyme not in enyzme_info.keys():
        raise ValueError('Unable to identify enzyme info for the specified enzyme.')

    return True


def validate_fidelity_data(fidelity_data: pd.DataFrame) -> bool:
    # validate the fidelity data
    if fidelity_data.shape[0] == 0 or fidelity_data.shape[1] == 0:
        raise ValueError('Empty fidelity data.')
    if fidelity_data.shape[0] != fidelity_data.shape[1]:
        raise ValueError('fidelity_data.shape[0] != fidelity_data.shape[1]')
    if 0 != sum(
            [make_rev_compliment(c) != i for idx, (c, i) in enumerate(zip(fidelity_data.columns, fidelity_data.index))]
    ):
        raise ValueError('Correct fusion site pairs are not at diagonal position!')

    return True


def generate_target_variants_from_muts(mutations_1dx: Optional[List[Any]], linked_mutations_1dx: Optional[List[Any]],
                                       wt_seq: str, incl_wt: bool) -> List:
    """
    """
    d: dict = {}
    linked_mut_positions = []
    if linked_mutations_1dx:
        for mut_set in linked_mutations_1dx:
            position_set = tuple(sorted([mut[1] for mut in mut_set]))
            mut_annotation_set = [f'{mut[0]}{mut[1]}{mut[2]}' for mut in mut_set]
            if position_set not in d.keys():
                d[position_set] = [mut_annotation_set]
            else:
                d[position_set].append(mut_annotation_set)
            if not incl_wt: continue
            d[position_set].append(['wt'])
        linked_mut_positions = [mut[1] for mut_set in linked_mutations_1dx for mut in mut_set]
    # print(d)
    if mutations_1dx:
        for mut in mutations_1dx:
            if mut['position'] in linked_mut_positions: continue
            d[mut['position']] = [[f"{wt_seq[mut['position'] - 1]}{mut['position']}{aa}"] for aa in mut['aa']]
            if not incl_wt: continue
            d[mut['position']].append(['wt'])
    all_var_notations = itertools.product(*[v for v in d.values()])
    var_notations_lst = []
    for var in all_var_notations:
        var_notations = '_'.join([mut for muts in var for mut in muts if mut != 'wt'])
        if var_notations != '': var_notations_lst.append(var_notations)
    return var_notations_lst


def compute_lib_complexity(mutations: Optional[List[Any]], linked_mutations: Optional[List[Any]]) -> int:
    """

    :param mutations: list of dictionary format
    :param linked_mutations: list of tuples format, must be in the identical indexing scheme as mutations
    :return:
    Example usage:
    mutations_1idx = include_linked_mutations_into_mutations(MUTATIONS, LINKED_MUTATIONS)
    compute_lib_complexity(mutations=mutations_1idx, linked_mutations=LINKED_MUTATIONS)
    """
    d: dict = {}
    linked_mut_positions = []
    if linked_mutations:
        for mut_set in linked_mutations:
            position_set = tuple(sorted([mut[1] for mut in mut_set]))
            if position_set not in d.keys():
                d.update({position_set: 1})
            else:
                d.update({position_set: d[position_set] + 1})
        linked_mut_positions = [mut[1] for mut_set in linked_mutations for mut in mut_set]
    # print(d)
    linked_mut_complexity = np.prod([v + 1 for v in d.values()])

    if mutations:
        independent_mut_complexity = int(math.prod(
            [len(mut['aa']) + 1 for mut in mutations if mut['position'] not in linked_mut_positions]
        ))
    else:
        independent_mut_complexity = 1
    return int(linked_mut_complexity * independent_mut_complexity)


@cache
def is_valid_fusion_site_set(fusion_sites: tuple) -> bool:
    # validate that each fusion site is unique, and no fusion sites are reverse compliment of each other
    fusion_sites_rev_comp = tuple(make_rev_compliment(i) for i in fusion_sites)
    if (len(fusion_sites) == len(set(fusion_sites))) and (
            len(fusion_sites_rev_comp + fusion_sites) == len(set(fusion_sites_rev_comp + fusion_sites))):
        return True
    return False


def is_dna(seq: str) -> bool:
    if re.fullmatch('[ATCG]+', seq):
        return True
    return False


def is_aa(seq: str) -> bool:
    if re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', seq):
        return True
    return False


def one_indexing_to_zero_indexing(mutations_1idx: list, linked_mutations_1idx: list) -> Tuple[list, list]:
    mutations_0_indexing = [{'position': mut['position'] - 1, 'aa': mut['aa']} for mut in mutations_1idx]
    linked_mutations_0idx = None
    if linked_mutations_1idx:
        linked_mutations_0idx = []
        for mut_set in linked_mutations_1idx:
            set_0idx = []
            for mut_1dx in mut_set:
                mut_0idx = (mut_1dx[0], mut_1dx[1] - 1) + mut_1dx[2:]
                set_0idx.append(mut_0idx)
            linked_mutations_0idx.append(tuple(set_0idx))

    return mutations_0_indexing, linked_mutations_0idx


def zero_indexing_to_one_indexing(mutations_0idx: Optional[List[Any]]) -> Optional[List[Dict[str, Any]]]:
    if mutations_0idx:
        mutations_1_indexing = [{'position': mut['position'] + 1, 'aa': mut['aa']} for mut in mutations_0idx]
    else:
        mutations_1_indexing = None

    return mutations_1_indexing


def frag_in_partition_too_short_or_too_long(tup: Union[tuple, list], min_aa_length: int, provider_max_dna_len: int,
                                            enzyme: str, enzyme_info_dic: Dict[str, dict]) -> bool:
    """
    Checks if the absolute difference between any two neighboring values in a list is greater than a given value N.

    Parameters:
    tup (tuple of int): A tup of integers.
    N (int): The threshold value to compare differences against.

    Returns:
    bool: True if any absolute difference between neighboring values in the list is greater than N, otherwise True.

    Example usage: numbers = (1, 4, 7, 10]) N = 3 print(frag_too_short(numbers, length_threshold))  # This will check
    if any neighboring pairs have a difference greater than 3
    """
    # In case of no cuts, just return PASS
    if len(tup) == 0:
        return True

    # Compute maximum aa fragment length
    fusion_site_len = enzyme_info_dic[enzyme]['fusion_site_length']
    recognition_site_len = len(enzyme_info_dic[enzyme]['recognition_site'])
    stuffer_len = 5
    max_aa_length = (provider_max_dna_len - fusion_site_len - (2 * (recognition_site_len + stuffer_len))) / 3

    for pos in range(1, len(tup)):
        frag_length = abs(tup[pos] - tup[pos - 1])
        if frag_length < min_aa_length or frag_length > max_aa_length:
            return True
    return False


def frag_in_partition_too_short(tup: Union[tuple, list], min_aa_length: int) -> bool:
    """
    Checks if the absolute difference between any two neighboring values in a list is greater than a given value N.

    Parameters:
    tup (tuple of int): A tup of integers.
    N (int): The threshold value to compare differences against.

    Returns:
    bool: True if any absolute difference between neighboring values in the list is greater than N, otherwise True.

    Example usage: numbers = (1, 4, 7, 10]) N = 3 print(frag_too_short(numbers, length_threshold))  # This will check
    if any neighboring pairs have a difference greater than 3
    """
    for pos in range(1, len(tup)):
        frag_length = abs(tup[pos] - tup[pos - 1])
        if frag_length < min_aa_length:
            return True
    return False


@cache
def make_rev_compliment(dna: str) -> str:
    s_rev = dna[::-1]
    compliment_nt_dic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    rev_comp = ''.join([compliment_nt_dic[s_rev[idx]] for idx in range(0, len(dna))])
    return rev_comp


def is_rev_compliment(s1: str, s2: str) -> bool:
    """
    Example usage:
    is_rev_compliment('ASGC','GATT')

    """
    compliment_nt_dic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for idx in range(0, len(s1)):
        if s1[idx] not in ['A', 'T', 'C', 'G'] or s2[-idx - 1] not in ['A', 'T', 'C', 'G']:
            raise ValueError('Detected non-canonical base!')

        if compliment_nt_dic[s1[idx]] != s2[-idx - 1]:
            return False

    return True


def remove_linked_mutations_from_all_mutations(all_mutations: list, linked_mutations: list) -> list:
    reformatted_linked_mutations = reformat_linked_mutations(linked_mutations)

    all_mutations_reformat = []
    for mut in all_mutations:
        for aa in mut['aa']:
            all_mutations_reformat.append({'position': mut['position'], 'aa': [aa]})
    all_mutations_reformat_filtered = []
    for mut in all_mutations_reformat:
        if mut not in reformatted_linked_mutations:
            all_mutations_reformat_filtered.append(mut)

    output = pd.DataFrame(all_mutations_reformat_filtered) \
        .groupby(['position']) \
        .agg({'aa': lambda x: list(itertools.chain.from_iterable(x))}) \
        .reset_index() \
        .to_dict(orient='records')

    return output


def reformat_linked_mutations(linked_mutations: list) -> list:
    """
    Example Input: [(('Q', 13, 'P'), ('P', 14, 'A')), (('R', 53, 'G'), ('V', 79, 'L'))]
    Example Output: [
    {'position': 13, 'aa': ['P']},
    {'position': 14, 'aa': ['A']},
    {'position': 53, 'aa': ['G']},
    {'position': 79, 'aa': ['L']}
    ]
    """
    reformat_muts = []
    if linked_mutations:
        for mut_set in linked_mutations:
            for mut in mut_set:
                reformat_muts.append({'position': mut[1], 'aa': [mut[2]]})

    return reformat_muts


def reformat_mutations(lst_of_tup_mutations: list) -> list:
    """
    Example Input: [('Q', 13, 'P'), ('A', 24, 'V'), ('G', 26, 'E'), ('R', 27, 'D'), ('T', 28, 'I'), ('T', 28, 'Q')]
    Example Output: [
    {'position': 13, 'aa': ['P']},
    {'position': 24, 'aa': ['V']},
    {'position': 26, 'aa': ['E']},
    {'position': 27, 'aa': ['D']},
    {'position': 28, 'aa': ['I', 'Q']}
    ]
    """
    reformat_muts = []
    seen_positions = []
    for mut in lst_of_tup_mutations:
        if mut[1] not in seen_positions:
            seen_positions.append(mut[1])
            reformat_muts.append({'position': mut[1], 'aa': [mut[2]]})
        else:
            sel_idx = [idx for idx, item in enumerate(reformat_muts) if item['position'] == mut[1]][0]
            reformat_muts[sel_idx].update({'aa': reformat_muts[sel_idx]['aa'] + [mut[2]]})
    reformat_muts = sorted(reformat_muts, key=lambda x: (x['position']))
    return reformat_muts


def include_linked_mutations_into_mutations(mutations: list, linked_mutations: list) -> list:
    """
    :param mutations: example: [
    {'position': 13, 'aa': ['P']},
    {'position': 24, 'aa': ['V']},
    {'position': 26, 'aa': ['E']},
    {'position': 27, 'aa': ['D']},
    {'position': 28, 'aa': ['I', 'Q']}
    ]
    :param linked_mutations: example: [
    {'position': 13, 'aa': ['P']},
    {'position': 14, 'aa': ['A']},
    {'position': 53, 'aa': ['G']},
    {'position': 79, 'aa': ['L']}
    ]
    :return: [
    {'position': 13, 'aa': ['P']},
    {'position': 14, 'aa': ['A']},
    {'position': 24, 'aa': ['V']},
    {'position': 26, 'aa': ['E']},
    {'position': 27, 'aa': ['D']},
    {'position': 28, 'aa': ['I', 'Q']},
    {'position': 53, 'aa': ['G']},
    {'position': 79, 'aa': ['L']}
    ]
    """
    reformatted_linked_mutations = reformat_linked_mutations(linked_mutations)
    seen_positions = []
    combined_mutations = []
    for mut in mutations + reformatted_linked_mutations:
        if mut['position'] not in seen_positions:
            seen_positions.append(mut['position'])
            combined_mutations.append({'position': mut['position'], 'aa': mut['aa']})
        else:
            sel_idx = [idx for idx, item in enumerate(combined_mutations) if item['position'] == mut['position']][0]
            combined_mutations[sel_idx].update({'aa': list(set(combined_mutations[sel_idx]['aa'] + mut['aa']))})
    combined_mutations = sorted(combined_mutations, key=lambda x: (x['position']))

    return combined_mutations


def generate_aa2codon_dict(codon_usage_tbl_path: str) -> dict:
    """"""
    codon_usage_tbl = pd.read_csv(codon_usage_tbl_path)
    codon_usage_dict_reformat = {}
    for idx, row in codon_usage_tbl.iterrows():
        if row['amino_acid'] not in codon_usage_dict_reformat:
            codon_usage_dict_reformat.update({
                row['amino_acid']: {
                    'codon': [row['codon']],
                    'relative_frequency': [row['relative_frequency']]
                }
            })
        else:
            codon_usage_dict_reformat[row['amino_acid']]['codon'].append(row['codon'])
            codon_usage_dict_reformat[row['amino_acid']]['relative_frequency'].append(row['relative_frequency'])

    return codon_usage_dict_reformat


def multi_well_plate_position_generator(row_range: Tuple[str, str], columns_range: Tuple[object, object]) -> list:
    """row will be capital alphabets
    example usage:
    multi_well_plate_position_generator(row_range=('A','D'), columns_range=(1,6))
    """

    # Initialize the character variable
    row = row_range[0]
    # Use a do-while loop to traverse and print the uppercase alphabets
    wells = []
    while True:
        for col in range(columns_range[0], columns_range[1] + 1):
            wells.append(''.join([row, str(col)]))
        row = chr(ord(row) + 1)
        if row > row_range[1]:
            break
    return wells


def find_best_codon_by_usage(codon_usage_tbl_path: str, aa: str, fix_base: dict) -> str:
    codon_usage_dict = generate_aa2codon_dict(codon_usage_tbl_path)
    max_fraction = -float('inf')
    sel_codon = "Not found"
    for codon, fraction in zip(codon_usage_dict[aa]['codon'], codon_usage_dict[aa]['relative_frequency']):
        passed = True
        for idx, base in fix_base.items():
            if base and codon[idx] != base:
                passed = False
        if passed and fraction > max_fraction:
            max_fraction = fraction
            sel_codon = codon
    return re.sub("U", "T", sel_codon)


def expand_python_codon_tables(name: str, taxid: int) -> None:
    """
    Expand the codon usage table of the external package: python-codon-tables
    https://pypi.org/project/python-codon-tables/
    """
    # name has to be in the format genus(one aphabet)_epithet, ex: e_coli, h_sapiens, c_griseus
    local_available_tables = [int(re.sub('^.*_', '', tbl)) for tbl in pct.available_codon_tables_names]
    if taxid not in local_available_tables:
        # Download it from online:
        # The data comes from http://www.kazusa.or.jp (they computed the codon usages from NCBI sequence data)
        table_dict = pct.get_codons_table(taxid)
        pct_dir = pct.__file__
        tbl_dir = path.join(path.dirname(path.dirname(pct_dir)), 'codon_usage_data', 'tables')
        tbl_name = "".join([name, '_', str(taxid), '.csv'])
        d = []
        for aa, codons in table_dict.items():
            for codon, relative_frequency in codons.items():
                d.append({'amino_acid': aa,
                          'codon': codon,
                          'relative_frequency': relative_frequency})
        pd.DataFrame(d).to_csv(path.join(tbl_dir, tbl_name), index=False)
        # remove blank line at the end of the file
        with open(path.join(tbl_dir, tbl_name), "r+", encoding="utf-8") as file:

            # Move the pointer (similar to a cursor in a text editor) to the end of the file
            file.seek(0, os.SEEK_END)

            # This code means the following code skips the very last character in the file -
            # i.e. in the case the last line is null we delete the last line
            # and the penultimate one
            pos = file.tell() - 1

            # Read each character in the file one at a time from the penultimate
            # character going backwards, searching for a newline character
            # If we find a new line, exit the search
            while pos > 0 and file.read(1) != "\n":
                pos -= 1
                file.seek(pos, os.SEEK_SET)

            # So long as we're not at the start of the file, delete all the characters ahead
            # of this position
            if pos > 0:
                file.seek(pos, os.SEEK_SET)
                file.truncate()

        orgs_tbl_path = path.join(tbl_dir, '..', 'organisms.csv')
        orgs_tbl = pd.read_csv(orgs_tbl_path)
        df = pd.concat([pd.DataFrame([[name, taxid]], columns=orgs_tbl.columns), orgs_tbl],
                       ignore_index=True).drop_duplicates()
        df.to_csv(orgs_tbl_path, index=False)


def find_first_stop_codon_idx(dna: str) -> Union[Tuple[int, int], None]:
    stop_codons = ['TAA', 'TAG', 'TGA']
    for pos in range(0, len(dna), 3):
        if dna[pos:pos + 3] in stop_codons:
            return pos, pos + 3
    return None


def prepare_0idx_mutations(mutations_1idx: List, linked_mutations_1idx: List) -> Tuple[List, List]:
    all_mutations_1idx = include_linked_mutations_into_mutations(mutations_1idx, linked_mutations_1idx)

    all_mutations_0idx, linked_mutations_0idx = one_indexing_to_zero_indexing(
        mutations_1idx=all_mutations_1idx, linked_mutations_1idx=linked_mutations_1idx
    )
    return all_mutations_0idx, linked_mutations_0idx


def breadth_first_product(*sequences: List) -> Generator:
    """Breadth First Search Cartesian Product
    Inspiration: https://stackoverflow.com/questions/42288203/generate-itertools-product-in-different-order
    """

    def partitions(n: int, k: int) -> Generator:
        for c in itertools.combinations(range(n + k - 1), k - 1):
            yield (b - a - 1 for a, b in zip((-1,) + c, c + (n + k - 1,)))

    max_position = [len(i) - 1 for i in sequences]
    for i in range(sum(max_position)):
        for positions in partitions(i, len(sequences)):
            try:
                yield tuple(map(lambda seq_, pos: seq_[pos], sequences, positions))
            except IndexError:
                continue
    yield tuple(map(lambda seq_, pos: seq_[pos], sequences, max_position))


def nearest_first_product(*sequences: List) -> Generator:
    """nearest_first_product.
    Inspiration: https://stackoverflow.com/questions/42288203/generate-itertools-product-in-different-order
    """
    start = (0,) * len(sequences)
    queue = [(0, start)]
    seen = {start}  # set(start)
    while queue:
        priority, indexes = heapq.heappop(queue)
        yield tuple(seq_[index] for seq_, index in zip(sequences, indexes))
        for i in range(len(sequences)):
            if indexes[i] < len(sequences[i]) - 1:
                lst = list(indexes)
                lst[i] += 1
                new_indexes = tuple(lst)
                if new_indexes not in seen:
                    new_priority = sum(index * index for index in new_indexes)
                    heapq.heappush(queue, (new_priority, new_indexes))
                    seen.add(new_indexes)


def depth_first_product(*sequences: List) -> Generator:
    for item in itertools.product(*sequences):
        yield item


if __name__ == '__main__':
    print('Some example usage of functions in utils.')
    fusion_sites_ = ['AG', 'CC', 'G', 'G']
    is_valid_fusion_site_set(fusion_sites_)

    mutations_ = [('Q', 13, 'P'), ('P', 14, 'A'), ('A', 24, 'V'), ('G', 26, 'E'), ('R', 27, 'D'), ('T', 28, 'I'),
                  ('T', 28, 'Q'), ('F', 29, 'M'), ('S', 30, 'D'), ('L', 45, 'N'), ('L', 47, 'Y'), ('R', 53, 'G'),
                  ('V', 79, 'L'), ('P', 100, 'K'), ('P', 100, 'R'), ('L', 109, 'F'), ('Y', 112, 'R')]
    reformat_mutations(lst_of_tup_mutations=mutations_)
