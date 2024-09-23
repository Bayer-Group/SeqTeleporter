"""Define function - Computes Fragment Length Unevenness"""


def compute_fragment_length_unevenness(s: str, partition: tuple) -> float:
    """
    Computes the unevenness in the lengths of fragments of a DNA sequence.

    This function calculates the unevenness in the lengths of fragments of a DNA sequence that is partitioned at
    specific indices. The unevenness is calculated as the ratio of the maximum fragment length to the minimum
    fragment length, subtracted by 1.

    Parameters:
    - s (str): The DNA sequence which is to be partitioned.
    - partition (tuple of int): Indices in the DNA sequence where it is partitioned. The start and end points
      of the sequence are automatically included as partition points.

    Returns:
    - float: The unevenness in the lengths of the fragments. The value is a ratio of the maximum fragment length
      to the minimum fragment length, subtracted by 1. A value of 0 indicates that all fragments are of equal length.

    Example:
    For a sequence "ATCGTACGTA" with partition points at (2, 5), the fragments are "AT", "CGT", "ACGTA". The
    unevenness is calculated as (5/2) - 1 = 1.5.
    """
    partition = (0,) + partition + (len(s),)
    fragment_sizes = [partition[idx] - partition[idx - 1] for idx in range(1, len(partition))]
    fragment_length_evenness = max(fragment_sizes) / min(fragment_sizes) - 1
    return fragment_length_evenness


if __name__ == "__main__":
    # Example usage:
    seq = 'SAEWTVEQDGMAIC'
    partition_ = (6, 10)
    compute_fragment_length_unevenness(s=seq, partition=partition_)
