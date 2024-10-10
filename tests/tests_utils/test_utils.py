from unittest import TestCase
import pandas as pd
from os import path, listdir
import re
from os.path import dirname, abspath
from python_codon_tables.python_codon_tables import _tables_dir

from seqteleporter.utils.utils import validate_fidelity_data, validate_enzyme_and_enzyme_info, compute_lib_complexity, \
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
        codon_usage_tbl_dir = _tables_dir
        host = 'c_griseus'
        codon_usage_tbl_path = ''
        for f in listdir(_tables_dir):
            if re.search(host, f):
                codon_usage_tbl_path = path.join(codon_usage_tbl_dir, f)
                break

        output = generate_aa2codon_dict(codon_usage_tbl_path)
        expected_output = {'*': {'codon': ['TAA', 'TAG', 'TGA'], 'relative_frequency': [0.26, 0.24, 0.5]},
                           'A': {'codon': ['GCA', 'GCC', 'GCG', 'GCT'], 'relative_frequency': [0.23, 0.37, 0.07, 0.32]},
                           'C': {'codon': ['TGC', 'TGT'], 'relative_frequency': [0.53, 0.47]},
                           'D': {'codon': ['GAC', 'GAT'], 'relative_frequency': [0.53, 0.47]},
                           'E': {'codon': ['GAA', 'GAG'], 'relative_frequency': [0.41, 0.59]},
                           'F': {'codon': ['TTC', 'TTT'], 'relative_frequency': [0.53, 0.47]},
                           'G': {'codon': ['GGA', 'GGC', 'GGG', 'GGT'], 'relative_frequency': [0.25, 0.34, 0.21, 0.2]},
                           'H': {'codon': ['CAC', 'CAT'], 'relative_frequency': [0.56, 0.44]},
                           'I': {'codon': ['ATA', 'ATC', 'ATT'], 'relative_frequency': [0.14, 0.51, 0.35]},
                           'K': {'codon': ['AAA', 'AAG'], 'relative_frequency': [0.39, 0.61]},
                           'L': {'codon': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'relative_frequency': [0.08, 0.19, 0.39, 0.13, 0.06, 0.14]},
                           'M': {'codon': ['ATG'], 'relative_frequency': [1.0]},
                           'N': {'codon': ['AAC', 'AAT'], 'relative_frequency': [0.55, 0.45]},
                           'P': {'codon': ['CCA', 'CCC', 'CCG', 'CCT'], 'relative_frequency': [0.29, 0.32, 0.08, 0.31]},
                           'Q': {'codon': ['CAA', 'CAG'], 'relative_frequency': [0.24, 0.76]},
                           'R': {'codon': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 'relative_frequency': [0.19, 0.19, 0.14, 0.18, 0.19, 0.11]},
                           'S': {'codon': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 'relative_frequency': [0.22, 0.15, 0.14, 0.22, 0.05, 0.22]},
                           'T': {'codon': ['ACA', 'ACC', 'ACG', 'ACT'], 'relative_frequency': [0.29, 0.37, 0.08, 0.26]},
                           'V': {'codon': ['GTA', 'GTC', 'GTG', 'GTT'], 'relative_frequency': [0.12, 0.24, 0.46, 0.18]},
                           'W': {'codon': ['TGG'], 'relative_frequency': [1.0]},
                           'Y': {'codon': ['TAC', 'TAT'], 'relative_frequency': [0.56, 0.44]}}
        self.assertEqual(expected_output, output)
