"""Define function - Get Cost of Ordering DNA"""

import math
from typing import Dict, Optional, List, Any


def find_cost(seq: str, partition: tuple, mutations_0idx: Optional[List[Any]],
              linked_mutations_0idx: Optional[List[Any]], enzyme: str, enzyme_info_dic: Dict[str, dict],
              cost_per_nt: float, provider_min_frag_len: int, provider_max_frag_len: int) -> float:

    """
    Calculates the total cost of synthesizing a DNA sequence (IDT eBlocks) with specific mutations using a given enzyme.

    The function considers the cost of synthesizing alternative fragments (mutations) within the sequence,
    incorporating the costs associated with fusion sites and recognition sites of the enzyme used for synthesis.
    It returns the total cost based on the number of nucleotides, including additional sequences required for
    the enzyme's action.

    Parameters:
    - seq (str): The DNA sequence for which the synthesis cost is to be calculated.
    - partition (tuple of int): Indices in the DNA sequence where it is partitioned for synthesis. These partitions
      are used to introduce mutations or for separate handling in the synthesis process.
    - mutations_0idx (list of dict): A list where each element is a dictionary representing a mutation. Each dictionary
      contains a 'position' key indicating the nucleotide position (0-indexed) of the mutation within the DNA sequence,
      and an 'aa' key with the amino acid(s) representing the mutation.
    - linked_mutations_0idx (list of tuple of int): A list where each element is a tuple of mutation positions that are
      linked. Each mutation position in a tuple is expected to occur simultaneously with the others.
    - enzyme (str): The name of the enzyme used for synthesis. This name is used to retrieve enzyme-specific information
      from `enzyme_info_dic`.
    - enzyme_info_dic (dict): A dictionary containing information about enzymes used in the synthesis process. Each key
      is an enzyme name, and its value is another dictionary with at least 'fusion_site_length' and 'recognition_site'
      keys, representing the length of the fusion site and the sequence of the recognition site, respectively.

    Returns:
    - float: The total cost of synthesizing the DNA sequence, including mutations, based on the cost per nucleotide
      and additional costs associated with the enzyme's fusion and recognition sites. If one of the fragment length is
      larger than max allowed provider specs, it returns NaN.

    Notes:
    - The cost per nucleotide is fixed at 0.06/bp (currency not user-adjustable).
    - The function calculates the cost of each fragment between partition points, considering mutations that fall within
      those fragments. It multiplies the number of alternative fragments (derived from the mutations) by the cost of
      synthesizing the fragment length in triplets (to account for codons), adding the costs for fusion and recognition
      sites. The total cost is the sum of costs for all fragments.
    """

    cost_sum = 0
    stuffer_len = 5
    partition_add_terminals = (0,) + partition + (len(seq),)
    linked_mutation_positions = []
    if linked_mutations_0idx:
        linked_mutation_positions = [[mut[1] for mut in mut_set] for mut_set in
                                     linked_mutations_0idx]  # example: [[13, 15],[35, 40]]
    for idx in range(0, len(partition_add_terminals) - 1):
        frag_variations = []
        frag_start, frag_end = partition_add_terminals[idx], partition_add_terminals[idx + 1]
        frag_len = frag_end - frag_start
        fusion_site_len = enzyme_info_dic[enzyme]['fusion_site_length']
        recognition_site_len = len(enzyme_info_dic[enzyme]['recognition_site'])
        frag_len_plus_recognition_site = (3 * frag_len) + fusion_site_len + (2 * (recognition_site_len + stuffer_len))
        if frag_len_plus_recognition_site < provider_min_frag_len:
            frag_len_plus_recognition_site = provider_min_frag_len
        if frag_len_plus_recognition_site > provider_max_frag_len and len(partition) > 0:  #
            # if one of the fragment length is larger than max allowed provider specs, return NaN (Not a Number)
            return float('nan')

        # print(frag_start,frag_end,frag_len,seq[frag_start:frag_end])
        if mutations_0idx:
            seen_linked_muts: list = []
            for mutation in mutations_0idx:
                is_linked_flag = False
                if frag_start <= mutation['position'] < frag_end:
                    for positions in linked_mutation_positions:
                        if mutation['position'] in positions:
                            is_linked_flag = True
                            if mutation['position'] not in seen_linked_muts:
                                # linked mutations have variation of 2: WT and all linked mutations presenting at the
                                # same time
                                frag_variations.append(2)
                                seen_linked_muts = seen_linked_muts + positions
                            break  # assume that one mutation location can only present in one linked mutation set
                    if not is_linked_flag:
                        frag_variations.append(len(mutation['aa']) + 1)

        n_alt_fragments = math.prod(frag_variations)
        sub_cost = n_alt_fragments * frag_len_plus_recognition_site * cost_per_nt
        cost_sum = cost_sum + sub_cost

    return round(cost_sum, 2)


def find_exact_cost(mutant_dna_fragments: list, price_per_base: float) -> float:
    total_cost = float(0)
    for frag in mutant_dna_fragments:
        frag_cost = len(frag.get('concatenated_dna', ""))*price_per_base
        total_cost += frag_cost

    return round(total_cost, 2)
