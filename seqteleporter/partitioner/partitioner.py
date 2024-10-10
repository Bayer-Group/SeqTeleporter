from itertools import combinations, chain, product
from typing import Union, List, Optional, Any, Dict

from seqteleporter.utils.utils import (include_linked_mutations_into_mutations, reformat_linked_mutations,
                                       breadth_first_product)


def find_cuttable_positions(s: str, mutations_0idx: Optional[List[Any]], linked_mutations_0idx: Optional[List[Any]],
                            min_aa_length: int, provider_max_dna_len: int, enzyme: str,
                            allowed_cut_positions_1idx: list, enzyme_info_dic: Dict[str, dict]) -> list[int]:
    # Compute maximum aa fragment length
    fusion_site_len = enzyme_info_dic[enzyme]['fusion_site_length']
    recognition_site_len = len(enzyme_info_dic[enzyme]['recognition_site'])
    stuffer_len = 5
    max_aa_len = (provider_max_dna_len - fusion_site_len - (2 * (recognition_site_len + stuffer_len))) / 3

    # List all cuttable_positions while respecting the minimum and maximum fragment length at both end
    n = len(s)
    cuttable_positions = [i for i in range(1, n) if
                          max_aa_len >= i >= min_aa_length and max_aa_len >= n - i >= min_aa_length]

    # Narrowing down cuttable_positions if user specifies allowed cut positions
    reforfmat_allowed_cut_positions: list = []
    for tup in allowed_cut_positions_1idx:
        if tup[1] != "end":
            reforfmat_allowed_cut_positions = reforfmat_allowed_cut_positions + list(range(tup[0] - 1, tup[1]))
        else:
            reforfmat_allowed_cut_positions = reforfmat_allowed_cut_positions + list(range(tup[0] - 1, n))
    if len(reforfmat_allowed_cut_positions) > 0:
        cuttable_positions = sorted(list(set(cuttable_positions).intersection(reforfmat_allowed_cut_positions)))

    # Narrowing down cuttable_positions: no cut adjacent to AA mutation site
    # ex: s="GGVVQPGRSL", mutation_site is at Proline, these two cuts are not allowed: "GGVVQ|PGRSL", "GGVVQP|GRSL"
    if mutations_0idx:
        mut_positions = sorted(list(set([mut['position'] for mut in mutations_0idx])))
    else:
        mut_positions = []
    for mut_position in mut_positions:
        if (min_aa_length <= mut_position < n - min_aa_length) and (mut_position in cuttable_positions) and (
                mut_position + 1 in cuttable_positions):
            cuttable_positions.remove(mut_position)
            cuttable_positions.remove(mut_position + 1)
        if mut_position == n - min_aa_length and (mut_position in cuttable_positions):
            cuttable_positions.remove(mut_position)
        # adapt to the case when consecutive mutation positions are specified (ex: position 12, 13)
        if mut_position + 1 in cuttable_positions:
            cuttable_positions.remove(mut_position + 1)
    # Narrowing down the cuttable_positions to allow linked mutations
    # Logic: if n mutations must appear together, then no cuts are allowed between the positions of these mutations.
    if linked_mutations_0idx:
        for linked_mutations in linked_mutations_0idx:
            positions = [mut[1] for mut in linked_mutations]
            start_pos = min(positions) + 1
            end_pos = max(positions)
            cuttable_positions = sorted(list(set(cuttable_positions) - set(range(start_pos, end_pos + 1))))

    return cuttable_positions


def find_even_cuts(string: str, regions: List) -> List:
    total_length = len(string)
    n = len(regions) + 1  # Number of parts will be n
    desired_length = total_length / n  # Desired length of each part

    cuts = []

    for i in range(len(regions)):
        start, end = regions[i]

        # Calculate the target length for this part
        target_cut_index = desired_length * (i + 1)

        # Find the best cut within the specified region
        if start <= target_cut_index <= end:
            cuts.append(int(round(target_cut_index,0)))
        elif target_cut_index < start:
            cuts.append(start)  # Cut at the beginning of the region
        else:
            cuts.append(end)  # Cut at the end of the region

    return cuts


def partitioner(s: str, cuttable_positions: list[int], number_of_cuts: int, mutations_0idx: list,
                linked_mutations_0idx: list, pre_distribute_mutations: bool, one_dist: bool) -> Union[list, str]:

    if number_of_cuts == 0:
        return [[tuple()]]

    # if there are not enough cuttable positions
    if len(cuttable_positions) < number_of_cuts:
        return (f'Too less cuttable positions. We need {number_of_cuts} cuts, '
                f'but we only have these cuttable positions: {cuttable_positions}')

    if not pre_distribute_mutations:
        partitions_list = [combinations(cuttable_positions, number_of_cuts)]
        return partitions_list

    # Find all the ways to distribute mutations and order them based on pseudo-cost (base count)
    mutation_distribution_dicts = distribute_mutations(s=s,
                                                       mutations_0idx=mutations_0idx,
                                                       linked_mutations_0idx=linked_mutations_0idx,
                                                       n_fragments=number_of_cuts + 1)
    # Generate cut ranges from each way to distribute mutations
    cut_ranges_list = generate_n_set_of_cut_ranges_from_a_list_of_mutation_distributions(
        mutation_distribution_dicts=mutation_distribution_dicts,
        top_n_set=len(mutation_distribution_dicts)
    )
    # number_of_cuts would be 0 when n < min_aa_length or when linked mutations flank entire sequence
    if len(cut_ranges_list) == 0:
        return (f'Can not pre distribute mutations. '
                f'Please try a different partition search strategy.')

    # Generate a list of partition generators,
    # each generator generates partitions for a way to distribut mutations
    partitions_list = []
    for cut_ranges in cut_ranges_list:
        sorted_cut_sites = sort_cut_sites_by_eveness(string=s, regions=cut_ranges)
        # make sure that all the cut sites are cuttable
        cut_positions = [list(set(lst).intersection(set(cuttable_positions))) for lst in sorted_cut_sites]
        # all of the cut ranges need to result in at least 1 cuttable positions
        if all([len(positions) > 0 for positions in cut_positions]):
            partitions_list.append(breadth_first_product(*cut_positions))
            if one_dist:  # if only one mutation distribution is used, then stop the loop
                break
    if len(partitions_list) == 0:
        return (f'Can not pre distribute mutations. '
                f'Please try a different partition search strategy.')

    return partitions_list


def sort_cut_sites_by_eveness(string: str, regions: List) -> List[List]:
    best_cuts = find_even_cuts(string, regions)
    sorted_cuttable_sites = sort_indices_by_distance(regions, best_cuts)
    return sorted_cuttable_sites


def sort_indices_by_distance(regions: List, best_cuts: List) -> List[List]:
    sorted_indices = []

    for i in range(len(regions)):
        region = list(range(regions[i][0],regions[i][1]+1))
        best_cut = best_cuts[i]

        # Calculate distances fromthe best cut index
        distances = [(abs(index - best_cut), index) for index in region]

        # Sort by distance (first element of tuple) and then by index if needed
        distances.sort(key=lambda x: (x[0], x[1]))

        # Append the sorted indices to the result
        sorted_indices.append([index for _, index in distances])

    return sorted_indices


def generate_cut_ranges_from_a_mutation_distribution(mutation_distribution: dict) -> list:
    allow_cut_ranges = []
    for idx, mutation_subset in enumerate(mutation_distribution['distribution']):
        if idx == len(mutation_distribution['distribution']) - 1:
            break
        cut_range_0idx = (mutation_subset[-1]['position'] + 2,
                          mutation_distribution['distribution'][idx + 1][0]['position'])
        # cut_range_1idx = tuple(i + 1 for i in cut_range_0idx)
        allow_cut_ranges.append(cut_range_0idx)

    return allow_cut_ranges


def count_bases_in_a_mutation_distribution(s: str, distributed_mutations_0idx_lists: List) -> int:
    base_counts = []
    for idx, muts in enumerate(distributed_mutations_0idx_lists):
        if idx == 0:
            frag_start = 0
            frag_end = (muts[-1]['position'] + distributed_mutations_0idx_lists[idx + 1][0]['position']) // 2
        elif idx == len(distributed_mutations_0idx_lists) - 1:
            frag_start = (muts[0]['position'] + distributed_mutations_0idx_lists[idx - 1][-1]['position']) // 2
            frag_end = len(s)
        else:
            frag_start = (muts[0]['position'] + distributed_mutations_0idx_lists[idx - 1][-1]['position']) // 2
            frag_end = (muts[-1]['position'] + distributed_mutations_0idx_lists[idx + 1][0]['position']) // 2
        frag_len = frag_end - frag_start
        fragment_count = 1
        for mut in muts:
            fragment_count *= len(mut['aa']) + 1
        base_count = fragment_count * frag_len
        base_counts.append(base_count)
    return sum(base_counts)


def distribute_mutations(s: str, mutations_0idx: list, linked_mutations_0idx: Optional[list],
                         n_fragments: int) -> List[dict]:
    all_mutations_0idx = include_linked_mutations_into_mutations(mutations=mutations_0idx,
                                                                 linked_mutations=linked_mutations_0idx)

    distribution_dicts = []
    begin, end = [0], [len(all_mutations_0idx)]
    if len(all_mutations_0idx) < n_fragments:
        raise ValueError("Unable to distribute mutations because there are less mutations than number of fragments.")

    mut_index_combinations = combinations(range(1, len(all_mutations_0idx)), n_fragments - 1)

    for distribution in mut_index_combinations:
        distributed = [all_mutations_0idx[sl] for sl in
                       map(slice, chain(begin, distribution), chain(distribution, end))]
        valid_distribution = True
        # Check for constraints
        for idx, mutation_subset in enumerate(distributed):
            if idx == len(distributed) - 1:
                break
            if linked_mutations_0idx:
                for muts in linked_mutations_0idx:
                    muts_lst = [muts]
                    reformatted_muts = reformat_linked_mutations(muts_lst)
                    if 0 < sum([mut in mutation_subset for mut in reformatted_muts]) < len(reformatted_muts):
                        # if only partial of a set of linked mutations locate within a group of mutation,
                        # the distribution is invalid
                        valid_distribution = False
                        break
            cut_range = distributed[idx + 1][0]['position'] - mutation_subset[-1]['position']
            # Constraint: In order to design fusion sites, we need at least 2 AAs between the last mutation position
            # of a distribution and the first mutation position of the next distribution.
            if cut_range <= 2:
                valid_distribution = False
                break

        if valid_distribution:
            base_count = count_bases_in_a_mutation_distribution(s, distributed)
            distribution_dicts.append({'distribution': distributed, 'base_count': base_count})

    sorted_indexes = sorted(range(len(distribution_dicts)), key=lambda i: distribution_dicts[i]['base_count'])

    return [distribution_dicts[i] for i in sorted_indexes]


def generate_n_set_of_cut_ranges_from_a_list_of_mutation_distributions(mutation_distribution_dicts: list,
                                                                       top_n_set: int) -> list[list]:
    sel_mutation_distribution_dicts = mutation_distribution_dicts[:top_n_set]
    if len(mutation_distribution_dicts) == 0:
        return []
    set_of_cut_ranges = []
    for mut_distribution_dict in sel_mutation_distribution_dicts:
        allow_cut_ranges = generate_cut_ranges_from_a_mutation_distribution(
            mutation_distribution=mut_distribution_dict
        )
        set_of_cut_ranges.append(allow_cut_ranges)
    return set_of_cut_ranges






