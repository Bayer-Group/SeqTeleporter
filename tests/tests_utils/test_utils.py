from unittest import TestCase
import pandas as pd
from os import path
from os.path import dirname, abspath

from proseqteleporter.utils.utils import validate_fidelity_data, validate_enzyme_and_enzyme_info, compute_lib_complexity, \
    validate_codon_table, generate_aa2codon_dict, is_dna, is_aa, pretty_fragments_expression


# example amino acid sequence is randomly generated
class TestPrettyFragmentsExpression(TestCase):

    def test_pretty_fragments_expression(self):
        fragments = [
            "NTHEIHWNRTEKDKP",
            "WRCSEN",
            "QVHGKFFPG",
            "TWEFYKCHCKTPHPL",
            "DWYSPYMMCSGDLKA",
            "YVAKYDPTHRTWWAV"
        ]
        fragment_with_fusion_sites = {
            "NTHEIHWNRTEKDKP": {
                "c_term_dna": "AA<ACCA|>",
                "middle_aa": "NTHEIHWNRTEKD"
            },
            "WRCSEN": {
                "n_term_dna": "<ACCA|>",
                "middle_aa": "WRCS",
                "c_term_dna": "GA<AAAC|>"
            },
            "QVHGKFFPG": {
                "n_term_dna": "<AAAC|>",
                "middle_aa": "QVHGKFFP",
                "c_term_dna": "GG<A|ACG>"
            },
            "TWEFYKCHCKTPHPL": {
                "n_term_dna": "<A|ACG>",
                "middle_aa": "WEFYKCHCKTPH",
                "c_term_dna": "CC<ACTC|>"
            },
            "DWYSPYMMCSGDLKA": {
                "n_term_dna": "<ACTC|>",
                "middle_aa": "DWYSPYMMCSGDL",
                "c_term_dna": "AA<AGCC|>"
            },
            "YVAKYDPTHRTWWAV": {
                "n_term_dna": "<AGCC|>",
                "middle_aa": "YVAKYDPTHRTWWAV"
            }
        }
        fusion_site_len = 4
        expr = pretty_fragments_expression(fragments, fragment_with_fusion_sites, fusion_site_len)
        expected_expr = \
            "\n[NTHEIHWNRTEKD]AA<ACCA|>" \
            "\n                 <ACCA|>[WRCS]GA<AAAC|>" \
            "\n                                <AAAC|>[QVHGKFFP]GG<A|ACG>" \
            "\n                                                   <A|ACG>[WEFYKCHCKTPH]CC<ACTC|>" \
            "\n                                                                          <ACTC|>[DWYSPYMMCSGDL]AA<AGCC|>" \
            "\n                                                                                                  <AGCC|>[YVAKYDPTHRTWWAV]"
        self.assertEqual(expected_expr, expr)


class TestIsDnaOrIsAA(TestCase):
    def test_is_dna_with_dna(self):
        s = 'ATCGAA'
        self.assertTrue(is_dna(s))

    def test_is_dna_without_dna(self):
        s = 'ADTCC'
        self.assertFalse(is_dna(s))

    def test_is_aa_with_aa(self):
        s = 'EFGHIKLMNPQRS'
        self.assertTrue(is_aa(s))

    def test_is_aa_without_aa(self):
        s = 'EFGHIKLMNZPQRS'
        self.assertFalse(is_aa(s))


class TestValidateCodonTable(TestCase):
    def test_validate_codon_table_valid_case(self):
        codon_table = {
            'A': ['GCA', 'GCC', 'GCG'],
            'C': ['TGT', 'TGC'],
            'D': ['GAT', 'GAC'],
            'E': ['GAA', 'GAG'],
            'F': ['TTT', 'TTC'],
            'G': ['GGT', 'GGC', 'GGA', 'GGG'],
            'H': ['CAT', 'CAC'],
            'I': ['ATT', 'ATC', 'ATA'],
            'K': ['AAA', 'AAG'],
            'L': ['CTT', 'CTC', 'CTG', 'TTA', 'TTG'],
            'M': ['ATG'],
            'N': ['AAT', 'AAC'],
            'P': ['CCT', 'CCC', 'CCA', 'CCG'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGT', 'CGC', 'CGG', 'AGA', 'AGG'],
            'S': ['TCT', 'TCC', 'TCG', 'AGT', 'AGC'],
            'T': ['ACT', 'ACC', 'ACA', 'ACG'],
            'V': ['GTT', 'GTC', 'GTA', 'GTG'],
            'W': ['TGG'],
            'Y': ['TAT', 'TAC']
        }
        self.assertTrue(validate_codon_table(codon_table))

    def test_validate_codon_table_invalid_case(self):
        codon_table = {}
        with self.assertRaises(ValueError) as context:
            validate_codon_table(codon_table)
        self.assertEqual('Provided codon table does not contain all 20 amino acids.', str(context.exception))


class TestValidateFidelityData(TestCase):
    def test_validate_fidelity_data_with_valid_input(self):
        fidelity_data = pd.DataFrame(
            {
                'GCCC': [259, 0, 0, 0],
                'AGTA': [0, 350, 0, 0],
                'GGGC': [0, 0, 259, 0],
                'TACT': [0, 0, 0, 350]
            }, index=['GGGC', 'TACT', 'GCCC', 'AGTA'])

        self.assertTrue(validate_fidelity_data(fidelity_data))

    def test_validate_fidelity_data_with_wrong_shape_input(self):
        fidelity_data = pd.DataFrame(
            {
                'GCCC': [259, 0, 0],
                'AGTA': [0, 350, 0],
                'GGGC': [0, 0, 259],
                'TACT': [0, 0, 0]
            }, index=['GGGC', 'TACT', 'GCCC'])

        with self.assertRaises(ValueError) as context:
            validate_fidelity_data(fidelity_data)
        self.assertEqual("fidelity_data.shape[0] != fidelity_data.shape[1]", str(context.exception))

    def test_validate_fidelity_data_with_wrong_index_input(self):
        fidelity_data = pd.DataFrame(
            {
                'GCCC': [259, 0, 0, 0],
                'AGTA': [0, 350, 0, 0],
                'GGGC': [0, 0, 259, 0],
                'TACT': [0, 0, 0, 350]
            }, index=['TACT', 'GGGC', 'GCCC', 'AGTA'])

        with self.assertRaises(ValueError) as context:
            validate_fidelity_data(fidelity_data)
        self.assertEqual(str(context.exception), "Correct fusion site pairs are not at diagonal position!")


class TestValidateEnzymeAndEnzymeInfo(TestCase):
    def setUp(self) -> None:
        self.enzyme_info = {
            'BsaI': {'fusion_site_length': 4,
                     'recognition_site': 'GGTCTCN'},
            'SapI': {'fusion_site_length': 3,
                     'recognition_site': 'GCTCTTCN'}
        }

    def test_validate_enzyme_and_enzyme_info_valid_case(self):
        enzyme = 'BsaI'
        self.assertTrue(validate_enzyme_and_enzyme_info(enzyme, self.enzyme_info))

    def test_validate_enzyme_and_enzyme_info_invalid_case(self):
        enzyme = 'BsaII'
        with self.assertRaises(ValueError) as context:
            validate_enzyme_and_enzyme_info(enzyme, self.enzyme_info)
        self.assertEqual('Unable to identify enzyme info for the specified enzyme.', str(context.exception))


class TestComputeLibComplexity(TestCase):
    def test_compute_lib_complexity_with_linked_mutations(self):
        mutations = [{'position': 5, 'aa': ['Q', 'S', 'T', 'A', 'D']},
                     {'position': 14, 'aa': ['T']},
                     {'position': 18, 'aa': ['T', 'M', 'A', 'P', 'K', 'S']},
                     {'position': 21, 'aa': ['M']},
                     {'position': 35, 'aa': ['G']}]
        linked_mutations = [(('G', 14, 'T'), ('C', 21, 'M'))]
        expected_output = 6 * 7 * 2 * 2
        self.assertEqual(expected_output, compute_lib_complexity(mutations, linked_mutations))


class TestGenerateAa2CodonDict(TestCase):
    def test_generate_aa2codon_dict(self):
        codon_usage_tbl_dir = path.join(dirname(dirname(dirname(abspath(__file__)))), 'proseqteleporter', 'data', 'codon_usage')
        host = 'c_griseus'
        output = generate_aa2codon_dict(codon_usage_tbl_dir, host)
        expected_output = {
            'F': {'triplet': ['TTT', 'TTC'], 'fraction': [0.47, 0.53]},
            'L': {'triplet': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                  'fraction': [0.06, 0.14, 0.13, 0.19, 0.08, 0.39]},
            'I': {'triplet': ['ATT', 'ATC', 'ATA'], 'fraction': [0.35, 0.51, 0.14]},
            'M': {'triplet': ['ATG'], 'fraction': [1.0]},
            'V': {'triplet': ['GTT', 'GTC', 'GTA', 'GTG'], 'fraction': [0.18, 0.24, 0.12, 0.46]},
            'S': {'triplet': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                  'fraction': [0.22, 0.22, 0.14, 0.05, 0.15, 0.22]},
            'P': {'triplet': ['CCT', 'CCC', 'CCA', 'CCG'], 'fraction': [0.31, 0.32, 0.29, 0.08]},
            'T': {'triplet': ['ACT', 'ACC', 'ACA', 'ACG'], 'fraction': [0.26, 0.37, 0.29, 0.08]},
            'A': {'triplet': ['GCT', 'GCC', 'GCA', 'GCG'], 'fraction': [0.32, 0.37, 0.23, 0.07]},
            'C': {'triplet': ['TGT', 'TGC'], 'fraction': [0.47, 0.53]},
            '*': {'triplet': ['TGA', 'TAA', 'TAG'], 'fraction': [0.5, 0.26, 0.24]},
            'W': {'triplet': ['TGG'], 'fraction': [1.0]},
            'R': {'triplet': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                  'fraction': [0.11, 0.18, 0.14, 0.19, 0.19, 0.19]},
            'G': {'triplet': ['GGT', 'GGC', 'GGA', 'GGG'], 'fraction': [0.2, 0.34, 0.25, 0.21]},
            'Y': {'triplet': ['TAT', 'TAC'], 'fraction': [0.44, 0.56]},
            'H': {'triplet': ['CAT', 'CAC'], 'fraction': [0.44, 0.56]},
            'Q': {'triplet': ['CAA', 'CAG'], 'fraction': [0.24, 0.76]},
            'N': {'triplet': ['AAT', 'AAC'], 'fraction': [0.45, 0.55]},
            'K': {'triplet': ['AAA', 'AAG'], 'fraction': [0.39, 0.61]},
            'D': {'triplet': ['GAT', 'GAC'], 'fraction': [0.47, 0.53]},
            'E': {'triplet': ['GAA', 'GAG'], 'fraction': [0.41, 0.59]}}
        self.assertEqual(expected_output, output)
