from unittest import TestCase
import math

from proseqteleporter.partition_property_finder.cost_finder import find_cost


class TestFindCost(TestCase):
    def setUp(self) -> None:
        self.partition = (10,)
        self.mutations_0idx = [{'position': 4, 'aa': ['Q', 'S', 'T', 'A', 'D']},
                          {'position': 13, 'aa': ['T', 'M', 'A', 'L']}, {'position': 17, 'aa': ['T']},
                          {'position': 20, 'aa': ['T', 'M']}, {'position': 34, 'aa': ['G']}]
        self.linked_mutations_0idx = [(('L', 17, 'P'), ('H', 34, 'A'))]
        self.enzyme = 'BsaI'
        self.enzyme_info_dic = {
            'BsaI': {'fusion_site_length': 4,
                     'recognition_site': 'GGTCTCN'},
            'SapI': {'fusion_site_length': 3,
                     'recognition_site': 'GCTCTTCN'}
        }
        self.cost_per_nt = 0.06
        self.provider_min_frag_len = 300
        self.provider_max_frag_len = 1500

    def test_find_cost_with_valid_seq_length(self):
        # example amino acid sequence is a part of published Human IgG1: https://rest.uniprot.org/uniprotkb/P0DOX5.fasta
        # first fragment < provider_min_frag_len (300 bp)
        # provider_max_frag_len > second fragment > provider_min_frag_len (300 bp)
        seq = "QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSL" \
              "RLSCAASGFTFSRYTIHWVRQA"

        expected_output = 752.4
        self.assertEqual(
            expected_output,
            find_cost(seq,
                      self.partition,
                      self.mutations_0idx,
                      self.linked_mutations_0idx,
                      self.enzyme,
                      self.enzyme_info_dic,
                      self.cost_per_nt,
                      self.provider_min_frag_len,
                      self.provider_max_frag_len)
        )

    def test_find_cost_exceeds_max_frag_len(self):
        # example amino acid sequence is a part of published Human IgG1: https://rest.uniprot.org/uniprotkb/P0DOX5.fasta
        seq = "QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSL" \
              "RLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHW" \
              "VRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQP" \
              "GRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRY" \
              "TIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAQVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQAVVQPGRSVVQ" \
              "PGRSVVQPGRSVVQPGRSVVQPGRSVVQPGRSVVQPGRSVVQPGRSVVQPGRSVVQPGRS"

        self.assertTrue(math.isnan(
            find_cost(seq,
                      self.partition,
                      self.mutations_0idx,
                      self.linked_mutations_0idx,
                      self.enzyme,
                      self.enzyme_info_dic,
                      self.cost_per_nt,
                      self.provider_min_frag_len,
                      self.provider_max_frag_len)
        ))

