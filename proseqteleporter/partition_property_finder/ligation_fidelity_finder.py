"""Define function - computes the ligation fidelity using NEB data"""
import pandas as pd
import numpy as np
from functools import cache
from typing import Dict
from os import environ, listdir
from os.path import join, dirname, abspath, basename
import static_frame as sf

from proseqteleporter.utils.utils import make_rev_compliment, validate_fidelity_data


# @cache
# def make_correct_lig_freq_dict(fidelity_data: sf.FrameHE) -> Dict[str, int]:
#     # if fidelity_data is not valid, validate_fidelity_data(fidelity_data) raises error, and code stops.
#     validate_fidelity_data(fidelity_data)
#     # transform to dict
#     sf_series = sf.Series(np.diag(fidelity_data))
#     correct_lig_freq_dict = {k: v for k, v in zip(fidelity_data.index, sf_series.values)}
#     # dict_ = pd.Series(np.diag(fidelity_data), index=[fidelity_data.index]).to_dict()
#     # correct_lig_freq_dict = {k[0]: v for k, v in dict_.items()}
#     return correct_lig_freq_dict


# @cache
# def compute_ligation_fidelity(all_fusion_sites_of_a_partition: tuple,
#                               fidelity_data: sf.FrameHE) -> float:
#     """
#     Computes the ligation fidelity for a partition of fusion sites. Fidelity is calculated as the
#     frequency of correct ligations divided by the total frequency of ligations for a given fusion site.
#     This involves calculating the fidelity of each junction within the partition, considering both the
#     original and reverse complement sites, and then multiplying these fidelities to get the total ligation
#     fidelity for the partition.
#
#     Parameters:
#     - all_fusion_sites_of_a_partition (list of str): Fusion sites within a partition, represented as strings.
#     - fidelity_data (pd.DataFrame): A DataFrame containing observed ligation frequencies between fusion sites.
#       Rows and columns represent fusion sites, and values indicate the frequency of ligation events.
#     - correct_lig_freq_dict (dict): A dictionary mapping fusion sites to their correct ligation frequencies,
#       as extracted from the diagonal of `fidelity_data`.
#
#     Returns:
#     - float: The total ligation fidelity for the partition, calculated as the product of individual
#       fidelities for each junction, taking into account both direct and reverse complement ligation events.
#
#     Notes:
#     - The function assumes that for any fusion site, the reverse complement site is also considered
#       in the fidelity calculation. If a site is its own reverse complement, special handling is applied
#       to avoid double counting.
#     """
#     validate_fidelity_data(fidelity_data)
#     # correct_lig_freq_dict = make_correct_lig_freq_dict(fidelity_data)
#     # transform to dict
#     sf_series = sf.Series(np.diag(fidelity_data))
#     correct_lig_freq_dict = {k: v for k, v in zip(fidelity_data.index, sf_series.values)}
#
#     all_fusion_sites_of_a_partition_set = set(all_fusion_sites_of_a_partition)  # remove duplicates
#     rev_comps = [make_rev_compliment(fs) for fs in all_fusion_sites_of_a_partition_set]
#     # remove duplicates that is present as reverse complement
#     sel_fusion_sites = list(set(list(all_fusion_sites_of_a_partition_set) + rev_comps))
#     sel_lig_freqs = fidelity_data.loc[sel_fusion_sites, sel_fusion_sites]
#     sel_lig_freqs_rowsums = sel_lig_freqs.sum()
#     total_lig_fidelity = 1
#     for fs in all_fusion_sites_of_a_partition_set:
#         rev_comp_fs = make_rev_compliment(fs)
#         if fs != rev_comp_fs:
#             fidelity_of_this_junction = \
#                 correct_lig_freq_dict[fs] * 2 / (
#                         sel_lig_freqs_rowsums[fs] + sel_lig_freqs_rowsums[rev_comp_fs])
#         else:
#             fidelity_of_this_junction = \
#                 correct_lig_freq_dict[fs] * 2 / (sel_lig_freqs_rowsums[fs] * 2 + correct_lig_freq_dict[fs] * 2)
#         total_lig_fidelity = total_lig_fidelity * fidelity_of_this_junction
#
#     return total_lig_fidelity



def compute_ligation_fidelity(all_fusion_sites_of_a_partition: tuple,
                              fidelity_data: np.ndarray,
                              fusion_site_indices: dict,
                              fusion_site_cols: dict,
                              correct_lig_freq_dict: dict) -> float:
    """
    Computes the ligation fidelity for a partition of fusion sites using a NumPy array for fidelity data.

    Parameters:
    - all_fusion_sites_of_a_partition (list of str): Fusion sites within a partition, represented as strings.
    - fidelity_data (np.ndarray): A 2D NumPy array containing observed ligation frequencies between fusion sites.
      The array is indexed by fusion sites using the provided fusion_site_indices.
    - fusion_site_indices (dict): A dictionary mapping fusion site names to their corresponding indices in the fidelity_data array.

    Returns:
    - float: The total ligation fidelity for the partition.
    """
    # Create a set of unique fusion sites and their reverse complements
    all_fusion_sites_of_a_partition_set = set(all_fusion_sites_of_a_partition)
    rev_comps = {make_rev_compliment(fs) for fs in all_fusion_sites_of_a_partition_set}
    sel_fusion_sites = list(all_fusion_sites_of_a_partition_set.union(rev_comps))

    # Get indices for the selected fusion sites
    sel_indices = [fusion_site_indices[fs] for fs in sel_fusion_sites]
    sel_indices_mapper = {fs: i for i, fs in enumerate(sel_fusion_sites)}
    sel_cols = [fusion_site_cols[fs] for fs in sel_fusion_sites]

    # Select ligation frequencies for the relevant fusion sites
    sel_lig_freqs = fidelity_data[np.ix_(sel_indices, sel_cols)]
    sel_lig_freqs_rowsums = sel_lig_freqs.sum(axis=1)

    total_lig_fidelity = 1.0

    for i, fs in enumerate(all_fusion_sites_of_a_partition_set):
        rev_comp_fs = make_rev_compliment(fs)
        fs_index_in_sel_lig_freqs = sel_indices_mapper[fs]
        if fs != rev_comp_fs:
            rev_comp_index_in_sel_lig_freqs = sel_indices_mapper[rev_comp_fs]
            fidelity_of_this_junction = (
                correct_lig_freq_dict[fs] * 2 /
                (sel_lig_freqs_rowsums[fs_index_in_sel_lig_freqs] + sel_lig_freqs_rowsums[rev_comp_index_in_sel_lig_freqs])
            )
        else:
            fidelity_of_this_junction = (
                correct_lig_freq_dict[fs] * 2 /
                (sel_lig_freqs_rowsums[fs_index_in_sel_lig_freqs] * 2 + correct_lig_freq_dict[fs] * 2)
            )
        total_lig_fidelity *= fidelity_of_this_junction

    return total_lig_fidelity


if __name__ == "__main__":
    # Example usage:
    fidelity_data_path = r'C:\Users\GOFKV\PycharmProjects\proseqteleporter\proseqteleporter\data\neb_fidelity_data\FileS01_T4_01h_25C.xlsx'
    fidelity_data_ = pd.read_excel(fidelity_data_path, index_col=0)
    all_fusion_sites_of_a_partition_ = ('AAGG', 'ACTC', 'AGGA', 'AGTG', 'ATCA')

    # Convert the DataFrame to a NumPy array
    fidelity_data = fidelity_data_.values

    # Create a mapping of fusion site names to their indices
    fusion_site_indices = {site: idx for idx, site in enumerate(fidelity_data_.index)}
    fusion_site_cols = {site: idx for idx, site in enumerate(fidelity_data_.columns)}
    res = compute_ligation_fidelity(all_fusion_sites_of_a_partition=all_fusion_sites_of_a_partition_,
                              fidelity_data=fidelity_data,
                              fusion_site_indices=fusion_site_indices,
                              fusion_site_cols=fusion_site_cols)


    # res = compute_ligation_fidelity(all_fusion_sites_of_a_partition=all_fusion_sites_of_a_partition_,
    #                                 fidelity_data=sf.FrameHE.from_pandas(fidelity_data_))
    print(res)