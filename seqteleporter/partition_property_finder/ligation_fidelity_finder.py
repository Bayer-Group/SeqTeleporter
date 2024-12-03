import numpy as np

from seqteleporter.utils.utils import make_rev_compliment


def compute_ligation_fidelity(all_fusion_sites_of_a_partition: tuple,
                              fidelity_data: np.ndarray,
                              fusion_site_indices: dict,
                              fusion_site_cols: dict,
                              correct_lig_freq_dict: dict) -> float:
    """
    Computes the ligation fidelity for a partition of fusion sites using a NumPy array for fidelity data.

    Parameters:
    - all_fusion_sites_of_a_partition (tuple of str): A tuple containing fusion site names within a partition.
    - fidelity_data (np.ndarray): A 2D NumPy array containing observed ligation frequencies between fusion sites.
      The array is indexed by fusion sites using the provided fusion_site_indices.
    - fusion_site_indices (dict): A dictionary mapping fusion site names to their corresponding indices in the fidelity_data array.
    - fusion_site_cols (dict): A dictionary mapping fusion site names to their respective column indices in the fidelity_data array.
    - correct_lig_freq_dict (dict): A dictionary containing the correct ligation frequencies for each fusion site.

    Returns:
    - float: The total ligation fidelity for the partition, calculated as the product of the fidelity for each junction.
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
