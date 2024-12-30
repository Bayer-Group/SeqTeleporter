from itertools import product
import numpy as np
import pandas as pd
import re
import math
from typing import Dict, List, Tuple, Optional, Any

from seqteleporter.utils.utils import (is_rev_compliment, frag_in_partition_too_short, is_valid_fusion_site_set,
                                       validate_fidelity_data, validate_enzyme_and_enzyme_info, validate_codon_table,
                                       depth_first_product, breadth_first_product, nearest_first_product)
from seqteleporter.partition_property_finder.ligation_fidelity_finder import compute_ligation_fidelity


class FusionSitesFinder:

    def assign_fusion_sites(
            self,
            s: str, mutations_0idx: Optional[List[Any]], partition: Tuple[int, ...], fidelity_data: pd.DataFrame,
            satisfaction_fidelity: float, enzyme: str, enzyme_info_dic: dict,
            fusion_sites_used_by_backbone: Tuple[str, ...], search_method: str, codon_table: Dict[str, List]
    ) -> Tuple[Optional[Tuple], float, List[List]]:
        """
        Assigns optimal fusion sites for each junction within a sequence based on experimental data, a fidelity
        threshold, and enzyme characteristics.

        This function evaluates potential fusion sites for sequence junctions created by a given partitioning. It
        leverages experimental data and enzyme information to identify fusion sites that meet or exceed a specified
        ligation fidelity threshold. The selection process generates candidate fusion sites for each junction,
        evaluates their ligation fidelity against the experimental data, and selects the set of fusion sites with the
        highest overall fidelity, or until it meets the satisfaction_fidelity criterion.

        Parameters:
            - s (str): The original sequence to be partitioned and fused at selected sites.
            - mutations_0idx (Optional[List[Any]]): A list of mutations indexed from 0, which may affect the fusion site selection.
            - partition (Tuple[int, ...]): Indices within the sequence `s` indicating where partitions (cuts) are made.
            - fidelity_data (pd.DataFrame): Experimental data used to calculate the ligation fidelity of potential fusion sites.
              This data should include information relevant to the efficiency of ligation at various sites within similar sequences.
            - satisfaction_fidelity (float): The fidelity threshold that must be met or exceeded for a fusion site to be selected.
            - enzyme (str): The name of the enzyme used in the partitioning process.
            - enzyme_info_dic (dict): A dictionary containing information about the enzymes used, including key characteristics like fusion site length.
            - fusion_sites_used_by_backbone (Tuple[str, ...]): A tuple of fusion sites that are reserved for backbone cloning and should not be considered.
            - search_method (str): The method used to search for optimal fusion sites. Options are 'DFS', 'BFS', or 'NFS'.
            - codon_table (Dict[str, List]): A mapping of codons to their corresponding amino acids, used for evaluating fusion site candidates.

        Returns:
            - Tuple[Optional[Tuple], float, List[List]]:
                1. sel_fusion_sites (tuple): The selected fusion sites for each junction that meet the satisfaction fidelity threshold.
                2. ligation_fidelity (float): The ligation fidelity of the selected fusion sites.
                3. sel_junction_dna_map_sliding_window (list of lists): Mapped information for each selected fusion site, detailing how they align with the sequence's underlying DNA and the experimental data.

        Raises:
            - ValueError: If any fragment resulting from the partitioning is shorter than 4 amino acids,
              as all fragments must meet this minimum length requirement to ensure viability.
            - ValueError: If an invalid search method is provided; must be one of 'NFS', 'BFS', or 'DFS'.
        """
        # In case of no cuts, just return default
        if len(partition) == 0:
            sel_fusion_sites: tuple = tuple()
            ligation_fidelity_of_sel_fusion_sites = float('nan')
            sel_junction_dna_map_sliding_window: list[list] = []
            return sel_fusion_sites, ligation_fidelity_of_sel_fusion_sites, sel_junction_dna_map_sliding_window

        validate_fidelity_data(fidelity_data)
        validate_enzyme_and_enzyme_info(enzyme, enzyme_info_dic)

        len_fusion_site = enzyme_info_dic[enzyme]['fusion_site_length']

        partition_with_terminals = list((0,) + partition + (len(s),))
        if frag_in_partition_too_short(tup=partition_with_terminals, min_aa_length=4):
            # validate all AA fragments length >= 4
            raise ValueError("At least One fragment too short!")

        candidate_fusion_sites_for_all_junctions = []
        junction_dna_map_sliding_window_all = []
        for cut in partition:
            junction_aa = s[cut - 2:cut + 2]
            # find_candidate_fusion_sites_for_a_junction
            unique_candidate_fusion_sites_for_this_junction, junction_dna_map_sliding_window = \
                self.find_candidate_fusion_sites_for_a_junction(junction_aa, len_fusion_site, codon_table)
            # refine_candidate_fusion_sites_for_a_cut
            refined_candidate_fusion_sites_for_this_junction, junction_dna_map_sliding_window = \
                self.refine_candidate_fusion_sites_for_a_cut(cut, mutations_0idx, junction_dna_map_sliding_window, codon_table)
            # remove fusions sites that are already reserved by user for backbone cloning
            refined_candidate_fusion_sites_for_this_junction = \
                [i for i in refined_candidate_fusion_sites_for_this_junction if i not in fusion_sites_used_by_backbone]
            # if any of the junction has no fitting fusion sites, this whole partition is invalid. return empty values.
            if len(refined_candidate_fusion_sites_for_this_junction) == 0:
                sel_fusion_sites: tuple = tuple()
                ligation_fidelity_of_sel_fusion_sites = float('nan')
                sel_junction_dna_map_sliding_window_all: list = []

                return sel_fusion_sites, ligation_fidelity_of_sel_fusion_sites, sel_junction_dna_map_sliding_window_all

            # collect refined fusion sites and junction_dna_map_sliding_window
            candidate_fusion_sites_for_all_junctions.append(refined_candidate_fusion_sites_for_this_junction)
            junction_dna_map_sliding_window_all.append(junction_dna_map_sliding_window)

        # select fusion sites with the best fidelity for each partition
        ligation_fidelity_of_sel_fusion_sites = -float('inf')
        sel_fusion_sites = tuple()
        loop = 0

        # DFS --------------------------------------------------------------------------------------------------------
        if search_method == 'DFS':
            candidate_fusion_sites_sets = depth_first_product(*candidate_fusion_sites_for_all_junctions)
        # BFS --------------------------------------------------------------------------------------------------------
        elif search_method == 'BFS':
            candidate_fusion_sites_sets = breadth_first_product(*candidate_fusion_sites_for_all_junctions)
        # NFS --------------------------------------------------------------------------------------------------------
        elif search_method == 'NFS':
            candidate_fusion_sites_sets = nearest_first_product(*candidate_fusion_sites_for_all_junctions)
        # ------------------------------------------------------------------------------------------------------------
        else:
            raise ValueError("Invalid search method. Search method must be one of these: 'NFS', 'BFS', 'DFS'.")
        # ------------------------------------------------------------------------------------------------------------
        # Create a mapping of fusion site names to their indices
        fusion_site_indices = {site: idx for idx, site in enumerate(fidelity_data.index)}
        fusion_site_cols = {site: idx for idx, site in enumerate(fidelity_data.columns)}
        # Extract correct ligation frequencies directly from the diagonal
        correct_lig_freq_dict = {site: fidelity_data.iloc[i, i] for site, i in fusion_site_indices.items()}
        for element in candidate_fusion_sites_sets:
            loop += 1
            all_fusion_sites = element + fusion_sites_used_by_backbone
            # remove cases where duplicated or rev-comp fusions sites are present in "element"
            if not is_valid_fusion_site_set(all_fusion_sites):
                continue
            ligation_fidelity = compute_ligation_fidelity(all_fusion_sites_of_a_partition=all_fusion_sites,
                                                          fidelity_data=fidelity_data.values,
                                                          fusion_site_indices=fusion_site_indices,
                                                          fusion_site_cols=fusion_site_cols,
                                                          correct_lig_freq_dict=correct_lig_freq_dict)
            if ligation_fidelity > ligation_fidelity_of_sel_fusion_sites:
                ligation_fidelity_of_sel_fusion_sites = ligation_fidelity
                sel_fusion_sites = element
            if ligation_fidelity_of_sel_fusion_sites >= satisfaction_fidelity:
                break

        # map to the junction_dna_map_sliding_window
        sel_junction_dna_map_sliding_window = []
        if sel_fusion_sites:
            for fusion_site, junction_dna_map_sliding_window in zip(sel_fusion_sites, junction_dna_map_sliding_window_all):
                sel_junction_dna_map_sliding_window.append(
                    [d for d in junction_dna_map_sliding_window if d['fusion_site'] == fusion_site]
                )

        return sel_fusion_sites, ligation_fidelity_of_sel_fusion_sites, sel_junction_dna_map_sliding_window

    def concat_sel_fusion_sites_to_fragments(
            self,
            fragments: List,
            fusion_sites: List,
            sel_junction_dna_map_fusion_sites: List,
            codon_usage_table_path: str
    ) -> dict:
        """
        Concatenates selected fusion sites to DNA fragments based on specific junctions and codon usage.

        This method processes a list of DNA fragments and fusion sites, determining how to
        concatenate them based on selected junctions mapped to the fusion sites. It generates
        an output dictionary that contains information about the order, C-terminal, and N-terminal
        DNA sequences for each fragment, alongside their corresponding middle amino acid sequences.

        Parameters:
        ----------
        fragments : List
            A list of DNA fragment identifiers (strings) to be processed.

        fusion_sites : List
            A list of fusion site sequences (strings) that will be concatenated to the fragments.

        sel_junction_dna_map_fusion_sites : List
            A list of selected junction DNA mappings corresponding to the fusion sites.

        codon_usage_table_path : str
            The file path to a codon usage table used for selecting the appropriate junction.

        Returns:
        -------
        dict
            A dictionary where each key is a fragment identifier and the value is another
            dictionary containing:
                - 'order': The order of the fragment in the concatenation process.
                - 'c_term_dna': The C-terminal DNA sequence after fusion.
                - 'n_term_dna': The N-terminal DNA sequence after fusion (if applicable).
                - 'middle_aa': The middle amino acid sequence after fusion.

        Notes:
        -----
        - If only one fragment is provided, the function returns a dictionary with that fragment's
          information without performing any fusion.
        - The fusion process involves determining how much of each fragment is retained and how
          much is replaced by the fusion site based on the selected junctions.
        """

        output_dict: dict = {}
        if len(fragments) == 1:
            output_dict[fragments[0]] = {'order': 0, 'middle_aa': fragments[0]}
            return output_dict

        for junction_idx, (fusion_site, junctions) in enumerate(zip(fusion_sites, sel_junction_dna_map_fusion_sites)):
            junction = self.select_junction_by_codon_usage(junctions, codon_usage_table_path)
            fs_start = junction['i']
            junction_dna = junction['junction_dna']
            cut_site = int(len(junction_dna) / 2)
            left_frag, right_frag = fragments[junction_idx], fragments[junction_idx + 1]

            if left_frag in output_dict:
                left_frag_edited = output_dict[left_frag]['middle_aa']
            else:
                left_frag_edited = left_frag

            min_i, max_i = len(junction_dna) / 2 - len(fusion_site), len(junction_dna) / 2
            left_frag_gives_away = math.ceil((len(fusion_site) - (fs_start - min_i)) / 3)
            fusion_site_expression = "".join([
                "<",
                junction_dna[fs_start:cut_site],
                "|",
                junction_dna[cut_site:fs_start + len(fusion_site)],
                ">"
            ])
            left_frag_c_term_dna = ''.join(
                [junction_dna[(2 - left_frag_gives_away) * 3:fs_start],
                 fusion_site_expression]
            )
            left_frag_middle_aa = left_frag_edited[:-left_frag_gives_away if left_frag_gives_away else None]
            right_frag_gives_away = math.ceil((len(fusion_site) - (max_i - fs_start)) / 3)
            right_frag_n_term_dna = ''.join([
                fusion_site_expression,
                junction_dna[
                    fs_start + len(fusion_site):-(2 - right_frag_gives_away) * 3 if 2 - right_frag_gives_away else None
                ]
            ])
            right_frag_middle_aa = right_frag[right_frag_gives_away:]
            if left_frag in output_dict:
                output_dict[left_frag].update({
                    'c_term_dna': left_frag_c_term_dna,
                    'middle_aa': left_frag_middle_aa
                })
                output_dict.update({
                    right_frag: {
                        'order': junction_idx + 1,
                        'n_term_dna': right_frag_n_term_dna,
                        'middle_aa': right_frag_middle_aa
                    }
                })
            else:
                output_dict.update({
                    left_frag: {
                        'order': junction_idx,
                        'c_term_dna': left_frag_c_term_dna,
                        'middle_aa': left_frag_middle_aa
                    },
                    right_frag: {
                        'order': junction_idx + 1,
                        'n_term_dna': right_frag_n_term_dna,
                        'middle_aa': right_frag_middle_aa
                    }
                })
        return output_dict

    @staticmethod
    def find_candidate_fusion_sites_for_a_junction(
            junction_aa: str,
            len_fusion_site: int,
            codon_table: Dict[str, List]
    ) -> Tuple[List, List]:
        """
        Identify potential DNA fusion sites for a given amino acid junction sequence.

        This function generates all possible DNA sequences corresponding to a junction amino acid sequence
        and then applies a sliding window approach to extract candidate fusion sites. The sliding window
        spans from (cut site - len_fusion_site) to (cut site + len_fusion_site), allowing for the
        identification of fusion sites based on the specified length.

        Parameters:
            junction_aa (str): The amino acid sequence at the junction for which candidate fusion sites are to be found.
            len_fusion_site (int): The length of the fusion site to be considered, which determines the size of the
                                   sliding window used to identify candidate fusion sites within the DNA sequence.
            codon_table (Dict[str, List]): A dictionary mapping single-letter amino acid codes to their corresponding
                                            codons (list of strings). This is used to translate the amino acid
                                            sequence into all possible DNA sequences.

        Returns:
            Tuple[List, List]: A tuple containing two elements:
                1. A list of unique DNA sequences representing potential fusion sites.
                2. A list of dictionaries, each representing a sliding window over the DNA sequence. Each dictionary
                   contains:
                    - 'junction_dna': The full DNA sequence corresponding to the junction.
                    - 'i': The start index of the fusion site within the DNA sequence.
                    - 'fusion_site': The DNA sequence of the fusion site itself.

        Notes:
            - The function removes palindromic fusion sites, as these are considered unusable in practical applications
              and can hinder processing speed.
        """
        validate_codon_table(codon_table)
        codon_lists = [codon_table[aa] for aa in junction_aa]
        unique_candidate_fusion_sites_for_this_junction: list = []
        junction_dna_map_sliding_window: list = []
        for junction_codons in product(*codon_lists):
            junction_dna = ''.join(junction_codons)
            cut_site = int(len(junction_dna) / 2)
            sliding_window_start = cut_site - len_fusion_site
            junction_dna_map_sliding_window = \
                junction_dna_map_sliding_window + \
                [{'junction_dna': junction_dna, 'i': idx, 'fusion_site': junction_dna[idx:idx + len_fusion_site]} for idx in
                 range(sliding_window_start, cut_site + 1)]
            sliding_windows = [junction_dna[idx:idx + len_fusion_site] for idx in range(sliding_window_start, cut_site + 1)]
            unique_candidate_fusion_sites_for_this_junction = sorted(
                list(set(unique_candidate_fusion_sites_for_this_junction + sliding_windows)))

        # Remove palindromic (double rotational symmetric) fusion sites as these sites are un-usable in realworld cases and
        # hinders process speed dramatically
        refined_candidate_fusion_sites_for_this_junction = []
        for fs in unique_candidate_fusion_sites_for_this_junction:
            if not is_rev_compliment(fs, fs):
                refined_candidate_fusion_sites_for_this_junction.append(fs)

        return refined_candidate_fusion_sites_for_this_junction, junction_dna_map_sliding_window

    @staticmethod
    def refine_candidate_fusion_sites_for_a_cut(cut: int, mutations_0idx: Optional[List[Any]],
                                                junction_dna_map_sliding_window: List,
                                                codon_table: Dict[str, List]) -> Tuple[List, List]:
        """
        Refines candidate fusion sites based on the provided cut position and mutations.

        This function filters candidate fusion sites at a given cut position by checking
        if any mutations overlap with the fusion site. It ensures that mutations do not
        lie within the fusion site and discards junctions where the mutations affect the
        candidate fusion sites.

        Parameters:
        - cut (int): The position of the cut in the DNA sequence (1-based index).
        - mutations_0idx (Optional[List[Any]]): A list of mutations, where each mutation
          is represented as a dictionary containing at least 'position' and 'aa' keys.
        - junction_dna_map_sliding_window (List): A list of dictionaries, each representing
          a candidate junction with keys 'junction_dna', 'i', and 'fusion_site'.
        - codon_table (Dict[str, List]): A dictionary mapping amino acid variations to
          their corresponding codons.

        Returns:
        - Tuple[List, List]: A tuple containing:
            1. A list of refined candidate fusion sites for the specified junction.
            2. The filtered junction DNA map sliding window, excluding discarded junctions.

        Raises:
        - ValueError: If no mutations are provided or if any mutation lies within the
          fusion site defined by the cut position.
        """
        # set default: dna codon of mutation position overlaps candidate fusion site at this junction
        discard_junction_dna_map_sliding_window = []
        if not mutations_0idx:
            raise ValueError('No desired mutations are given')
        else:
            for mut in mutations_0idx:
                # validate the cut position
                if mut['position'] == cut - 1 or mut['position'] == cut:
                    raise ValueError(f"Invalid cut position!\n "
                                     f"Cut position: {cut} results in mutation {mut} lying completely within fusion site")
                # dna codon of mutation position overlaps candidate fusion site at this junction
                else:
                    # mutation position is the second amino acid left to cut site
                    if mut['position'] == cut - 2:
                        for d in junction_dna_map_sliding_window:
                            junction_dna = d['junction_dna']
                            idx = d['i']
                            # dna codon of mutation position overlaps candidate fusion site
                            if idx < 3:
                                overlapped_bases = junction_dna[idx:3]
                                desired_variations = mut['aa']
                                for desired_variation in desired_variations:
                                    sub_codons_lst = [codon[-len(overlapped_bases):]
                                                      for codon in codon_table[desired_variation]]
                                    if overlapped_bases not in sub_codons_lst:
                                        discard_junction_dna_map_sliding_window.append(d)
                                        break
                    # mutation position is the second amino acid right to cut site
                    if mut['position'] == cut + 1:
                        for d in junction_dna_map_sliding_window:
                            idx = d['i']
                            junction_dna = d['junction_dna']
                            fs_length = len(d['fusion_site'])
                            # dna codon of mutation position overlaps candidate fusion site
                            if idx + fs_length > len(junction_dna) - 3:
                                overlapped_bases = junction_dna[-3:idx + fs_length]
                                desired_variations = mut['aa']
                                for desired_variation in desired_variations:
                                    sub_codons_lst = [codon[:len(overlapped_bases)]
                                                      for codon in codon_table[desired_variation]]
                                    if overlapped_bases not in sub_codons_lst:
                                        discard_junction_dna_map_sliding_window.append(d)
                                        break

        junction_dna_map_sliding_window = [
            d for d in junction_dna_map_sliding_window if d not in discard_junction_dna_map_sliding_window
        ]
        candidate_fusion_sites_for_this_junction = sorted(
            list(
                set(
                    [d['fusion_site'] for d in junction_dna_map_sliding_window]
                )
            )
        )
        return candidate_fusion_sites_for_this_junction, junction_dna_map_sliding_window

    @staticmethod
    def select_junction_by_codon_usage(junctions: List, codon_usage_table_path: str) -> dict:
        """Select junction DNA based on codon usage table"""

        # validate inputs
        if len(junctions) == 0:
            raise ValueError("No junctions are provided!")
        codon_usage_dict = pd.read_csv(codon_usage_table_path, index_col='codon').to_dict(orient='index')
        codon_usage_dict = {re.sub('U','T',key):val for key, val in codon_usage_dict.items()}
        max_score = -float('inf')
        sel_juction = {}
        for junction in junctions:
            junction_dna = junction['junction_dna']
            codon_usage_score = np.prod(
                [codon_usage_dict[junction_dna[i:i + 3]]['relative_frequency'] for i in range(0, len(junction_dna), 3)])
            if codon_usage_score > max_score:
                sel_juction = junction
        return sel_juction

