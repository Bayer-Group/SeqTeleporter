from unittest import TestCase

from proseqteleporter.partition_property_finder.fragment_length_uneveness_finder import compute_fragment_length_unevenness


class TestComputeFragmentLengthUnevenness(TestCase):

    def test_find_cost_with_valid_seq_length(self):
        # example amino acid sequence is a part of published Human IgG1: https://rest.uniprot.org/uniprotkb/P0DOX5.fasta
        seq = "QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSL" \
              "RLSCAASGFTFSRYTIHWVRQA"
        partition = (10, 33)

        expected_output = 87 / 10 - 1

        self.assertEqual(
            expected_output,
            compute_fragment_length_unevenness(s=seq, partition=partition)
        )
