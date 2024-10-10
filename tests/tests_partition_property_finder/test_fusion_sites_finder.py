import re
from unittest import TestCase
import pandas as pd
from os.path import dirname, abspath, join
from os import mkdir, listdir, path
import math
import shutil
from python_codon_tables.python_codon_tables import _tables_dir as codon_usage_table_dir

from proseqteleporter.partition_property_finder.fusion_sites_finder import \
    breadth_first_product, nearest_first_product, find_candidate_fusion_sites_for_a_junction, \
    refine_candidate_fusion_sites_for_a_cut, assign_fusion_sites, select_junction_by_codon_usage, \
    concat_sel_fusion_sites_to_fragments

from proseqteleporter.config import CODON_TABLE, ENZYME_INFO

FIDELITY_DATA_PATH = path.join(
    path.dirname(path.dirname(dirname(abspath(__file__)))),
    'proseqteleporter', 'data', 'neb_fidelity_data', 'FileS01_T4_01h_25C.xlsx'
)

host = 'c_griseus'
codon_usage_table_path = ""
for f in listdir(codon_usage_table_dir):
    if re.search(host, f):
        codon_usage_table_path = join(codon_usage_table_dir, f)
        break


class TestBreadthFirstProduct(TestCase):
    def test_breadth_first_product_normal_case(self):
        input_list = [list(range(1, 5)), list(range(1, 4)), list(range(1, 5))]
        expected_outputs = [(1, 1, 1), (1, 1, 2), (1, 2, 1), (2, 1, 1), (1, 1, 3), (1, 2, 2), (1, 3, 1), (2, 1, 2),
                            (2, 2, 1), (3, 1, 1), (1, 1, 4), (1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 2, 2), (2, 3, 1),
                            (3, 1, 2), (3, 2, 1), (4, 1, 1), (1, 2, 4), (1, 3, 3), (2, 1, 4), (2, 2, 3), (2, 3, 2),
                            (3, 1, 3), (3, 2, 2), (3, 3, 1), (4, 1, 2), (4, 2, 1), (1, 3, 4), (2, 2, 4), (2, 3, 3),
                            (3, 1, 4), (3, 2, 3), (3, 3, 2), (4, 1, 3), (4, 2, 2), (4, 3, 1), (2, 3, 4), (3, 2, 4),
                            (3, 3, 3), (4, 1, 4), (4, 2, 3), (4, 3, 2), (3, 3, 4), (4, 2, 4), (4, 3, 3), (4, 3, 4)]
        self.assertEqual(expected_outputs, [out for out in breadth_first_product(*input_list)])


class TestNearestFirstProduct(TestCase):
    def test_nearest_first_product_normal_case(self):
        input_list = [list(range(1, 5)), list(range(1, 4)), list(range(1, 5))]
        expected_outputs = [(1, 1, 1), (1, 1, 2), (1, 2, 1), (2, 1, 1), (1, 2, 2), (2, 1, 2), (2, 2, 1), (2, 2, 2),
                            (1, 1, 3), (1, 3, 1), (3, 1, 1), (1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2),
                            (3, 2, 1), (2, 2, 3), (2, 3, 2), (3, 2, 2), (1, 3, 3), (3, 1, 3), (3, 3, 1), (1, 1, 4),
                            (2, 3, 3), (3, 2, 3), (3, 3, 2), (4, 1, 1), (1, 2, 4), (2, 1, 4), (4, 1, 2), (4, 2, 1),
                            (2, 2, 4), (4, 2, 2), (3, 3, 3), (1, 3, 4), (3, 1, 4), (4, 1, 3), (4, 3, 1), (2, 3, 4),
                            (3, 2, 4), (4, 2, 3), (4, 3, 2), (3, 3, 4), (4, 3, 3), (4, 1, 4), (4, 2, 4), (4, 3, 4)]
        self.assertEqual(expected_outputs, [out for out in nearest_first_product(*input_list)])


class TestFindCandidateFusionSitesForAJunction(TestCase):
    def setUp(self) -> None:
        self.expected_output_junction_dna_map_sliding_window = [
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 5, 'fusion_site': 'GAAA'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 6, 'fusion_site': 'AAAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 5, 'fusion_site': 'GAAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 6, 'fusion_site': 'AAAA'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 5, 'fusion_site': 'GAAG'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 6, 'fusion_site': 'AAGA'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 5, 'fusion_site': 'GAAG'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 6, 'fusion_site': 'AAGA'},
        ]

    def test_find_candidate_fusion_sites_for_a_junction_simple(self):
        junction_aa = 'WWKK'
        len_fusion_site = 4
        unique_candidate_fusion_sites_for_this_junction, junction_dna_map_sliding_window = \
            find_candidate_fusion_sites_for_a_junction(junction_aa=junction_aa, len_fusion_site=len_fusion_site,
                                                       codon_table=CODON_TABLE)
        expected_output_unique_candidate_fusion_sites = {'GTGG', 'TGGA', 'GGAA', 'GAAA', 'AAAA', 'GAAG', 'AAGA'}
        expected_output = (expected_output_unique_candidate_fusion_sites,
                           self.expected_output_junction_dna_map_sliding_window)
        self.assertEqual(
            expected_output,
            (set(unique_candidate_fusion_sites_for_this_junction), junction_dna_map_sliding_window)
        )

    def test_find_candidate_fusion_sites_for_a_junction_with_palindromic_fusion_sites(self):
        junction_aa = 'MDHM'
        len_fusion_site = 4
        unique_candidate_fusion_sites_for_this_junction, junction_dna_map_sliding_window = \
            find_candidate_fusion_sites_for_a_junction(junction_aa=junction_aa, len_fusion_site=len_fusion_site,
                                                       codon_table=CODON_TABLE)
        expected_output_unique_candidate_fusion_sites = {'GGAT', 'ATCA', 'TCAT', 'CATA', 'TCAC', 'CCAT', 'GGAC',
                                                         'GACC', 'ACCA', 'CCAC', 'CACA'}
        expected_output_junction_dna_map_sliding_window = [
            {'junction_dna': 'ATGGATCATATG', 'i': 2, 'fusion_site': 'GGAT'},
            {'junction_dna': 'ATGGATCATATG', 'i': 3, 'fusion_site': 'GATC'},
            {'junction_dna': 'ATGGATCATATG', 'i': 4, 'fusion_site': 'ATCA'},
            {'junction_dna': 'ATGGATCATATG', 'i': 5, 'fusion_site': 'TCAT'},
            {'junction_dna': 'ATGGATCATATG', 'i': 6, 'fusion_site': 'CATA'},
            {'junction_dna': 'ATGGATCACATG', 'i': 2, 'fusion_site': 'GGAT'},
            {'junction_dna': 'ATGGATCACATG', 'i': 3, 'fusion_site': 'GATC'},
            {'junction_dna': 'ATGGATCACATG', 'i': 4, 'fusion_site': 'ATCA'},
            {'junction_dna': 'ATGGATCACATG', 'i': 5, 'fusion_site': 'TCAC'},
            {'junction_dna': 'ATGGATCACATG', 'i': 6, 'fusion_site': 'CACA'},
            {'junction_dna': 'ATGGACCATATG', 'i': 2, 'fusion_site': 'GGAC'},
            {'junction_dna': 'ATGGACCATATG', 'i': 3, 'fusion_site': 'GACC'},
            {'junction_dna': 'ATGGACCATATG', 'i': 4, 'fusion_site': 'ACCA'},
            {'junction_dna': 'ATGGACCATATG', 'i': 5, 'fusion_site': 'CCAT'},
            {'junction_dna': 'ATGGACCATATG', 'i': 6, 'fusion_site': 'CATA'},
            {'junction_dna': 'ATGGACCACATG', 'i': 2, 'fusion_site': 'GGAC'},
            {'junction_dna': 'ATGGACCACATG', 'i': 3, 'fusion_site': 'GACC'},
            {'junction_dna': 'ATGGACCACATG', 'i': 4, 'fusion_site': 'ACCA'},
            {'junction_dna': 'ATGGACCACATG', 'i': 5, 'fusion_site': 'CCAC'},
            {'junction_dna': 'ATGGACCACATG', 'i': 6, 'fusion_site': 'CACA'},
        ]
        expected_output = (expected_output_unique_candidate_fusion_sites,
                           expected_output_junction_dna_map_sliding_window)
        self.assertEqual(
            expected_output,
            (set(unique_candidate_fusion_sites_for_this_junction), junction_dna_map_sliding_window)
        )


class TestRefineCandidateFusionSitesForAPartition(TestCase):
    def setUp(self) -> None:
        self.cuts = [30, 27]  # test two overlap cases: mut['position'] == cut - 2, and mut['position'] == cut + 1
        self.junction_dna_map_sliding_window = [
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 5, 'fusion_site': 'GAAA'},
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 6, 'fusion_site': 'AAAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 5, 'fusion_site': 'GAAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 6, 'fusion_site': 'AAAA'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 5, 'fusion_site': 'GAAG'},
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 6, 'fusion_site': 'AAGA'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 2, 'fusion_site': 'GTGG'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 3, 'fusion_site': 'TGGA'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 4, 'fusion_site': 'GGAA'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 5, 'fusion_site': 'GAAG'},
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 6, 'fusion_site': 'AAGA'},
        ]
        self.candidate_fusion_sites_for_this_junction = ['GTGG', 'TGGA', 'GGAA', 'GAAA', 'AAAA', 'GAAG', 'AAGA']

    def synthesize_expected_and_actual_outputs(self, rmv_idxes, mutations_0idx):
        expected_outputs, actual_outputs = [], []
        for cut, rmv_idx in zip(self.cuts, rmv_idxes):
            expected_refined_junction_dna_map_sliding_window = [
                d for d in self.junction_dna_map_sliding_window if d['i'] != rmv_idx
            ]
            expected_refined_unique_candidate_fusion_sites = {
                d['fusion_site'] for d in expected_refined_junction_dna_map_sliding_window
            }
            expected_outputs.append(
                (expected_refined_unique_candidate_fusion_sites, expected_refined_junction_dna_map_sliding_window)
            )

            candidate_fusion_sites_for_this_junction_out, junction_dna_map_sliding_window_out = \
                refine_candidate_fusion_sites_for_a_cut(
                    cut,
                    mutations_0idx,
                    self.junction_dna_map_sliding_window,
                    CODON_TABLE
                )
            actual_outputs.append(
                (set(candidate_fusion_sites_for_this_junction_out), junction_dna_map_sliding_window_out)
            )
        return expected_outputs, actual_outputs

    def test_refine_candidate_fusion_sites_for_a_cut_with_fs_overlap_mut_and_fs_not_pass(self):
        mutations_0idx = [{'position': 28, 'aa': ['F', 'M']}]  # the fusion site fits for M, but not for F
        remove_idxes = [2, 6]  # test two overlap cases: mut['position'] == cut - 2, and mut['position'] == cut + 1
        expected_outputs, actual_outputs = self.synthesize_expected_and_actual_outputs(remove_idxes, mutations_0idx)
        for expected_output, actual_output in zip(expected_outputs, actual_outputs):
            self.assertEqual(expected_output, actual_output)

    def test_refine_candidate_fusion_sites_for_a_cut_with_fs_overlap_mut_and_fs_pass(self):
        mutations_0idx = [{'position': 28, 'aa': ['M', 'K', 'R']}]  # the fusion site fits for M, K, R
        remove_idxes = [None, None]  # no entries need to be removed for these two overlap cases
        expected_outputs, actual_outputs = self.synthesize_expected_and_actual_outputs(remove_idxes, mutations_0idx)
        for expected_output, actual_output in zip(expected_outputs, actual_outputs):
            self.assertEqual(expected_output, actual_output)

    def test_refine_candidate_fusion_sites_for_a_cut_without_fs_overlap_mut(self):
        mutations_0idx = [{'position': 10, 'aa': ['M', 'K']}]  # the fusion sites don't overlap mutation site
        remove_idxes = [None, None]  # no entries need to be removed for these two overlap cases
        expected_outputs, actual_outputs = self.synthesize_expected_and_actual_outputs(remove_idxes, mutations_0idx)
        for expected_output, actual_output in zip(expected_outputs, actual_outputs):
            self.assertEqual(expected_output, actual_output)


class TestAssignFusionSites(TestCase):
    def setUp(self) -> None:
        self.output_dir = path.join(dirname(abspath(__file__)), 'output_secret')

        self.inputs = dict(
            s='QVQLVQSWWKKVVQP',
            partition=(9,),
            mutations_0idx=[{'position': 7, 'aa': ['K']}, {'position': 10, 'aa': ['P']}],
            fidelity_data=pd.read_excel(FIDELITY_DATA_PATH, index_col=0),
            satisfaction_fidelity=0.95,
            enzyme='BsaI',
            enzyme_info_dic=ENZYME_INFO,
            fusion_sites_used_by_backbone=('CCAT', 'TGTC', 'TAAT'),
            search_method="BFS",
            codon_table=CODON_TABLE
        )
        if not path.exists(self.output_dir):
            mkdir(self.output_dir)

    def test_assign_fusion_sites_normal_cases(self):

        expected_sel_fusion_sites = ('GAAA',)
        expected_ligation_fidelity_of_sel_fusion_sites = 0.996
        expected_sel_junction_dna_map_sliding_window = [[
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 5, 'fusion_site': 'GAAA'},
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 5, 'fusion_site': 'GAAA'}
        ]]

        sel_fusion_sites_, ligation_fidelity_of_sel_fusion_sites_, sel_junction_dna_map_sliding_window_ = \
            assign_fusion_sites(**self.inputs)

        self.assertEqual(expected_sel_fusion_sites, sel_fusion_sites_)
        self.assertEqual(expected_ligation_fidelity_of_sel_fusion_sites,
                         round(ligation_fidelity_of_sel_fusion_sites_, 3))
        self.assertEqual(expected_sel_junction_dna_map_sliding_window, sel_junction_dna_map_sliding_window_)
        shutil.rmtree(self.output_dir)

    def test_assign_fusion_sites_edge_cases(self):

        edge_cases = [('s', ''), ('partition', tuple()), ('mutations_0idx', []), ('fidelity_data', pd.DataFrame({})),
                      ('enzyme', ''), ('enzyme_info_dic', {}), ('search_method', ''), ('codon_table', {})]

        for param, value in edge_cases:
            inputs = self.inputs.copy()
            inputs.update({param: value})
            expected_exceptions = {
                'codon_table': 'Provided codon table does not contain all 20 amino acids.',
                'enzyme': 'Unable to identify enzyme info for the specified enzyme.',
                'enzyme_info_dict': 'Unable to identify enzyme info for the specified enzyme.',
                'fidelity_data': 'Empty fidelity data.',
                'mutations_0idx': 'No desired mutations are given'
            }
            expected_outputs = {
                's': (tuple(), float('nan'), []),
                'partition': (tuple(), float('nan'), [])
            }

            if param in expected_outputs.keys():
                sel_fusion_sites, ligation_fidelity_of_sel_fusion_sites, sel_junction_dna_map_sliding_window = \
                    assign_fusion_sites(**inputs)
                self.assertEqual(expected_outputs[param][0], sel_fusion_sites)
                ligation_fidelity_of_sel_fusion_sites = round(ligation_fidelity_of_sel_fusion_sites, 3)
                if math.isnan(expected_outputs[param][1]):
                    self.assertTrue(math.isnan(ligation_fidelity_of_sel_fusion_sites))
                else:
                    self.assertEqual(expected_outputs[param][1], ligation_fidelity_of_sel_fusion_sites)
                self.assertEqual(expected_outputs[param][2], sel_junction_dna_map_sliding_window)
            if param in expected_exceptions.keys():
                with self.assertRaises(ValueError) as context:
                    assign_fusion_sites(**inputs)
                self.assertEqual(expected_exceptions[param], str(context.exception))
        shutil.rmtree(self.output_dir)

    def test_assign_fusion_sites_invalid_fidelity_data(self):
        inputs = self.inputs.copy()
        invalid_fidelity_data = pd.DataFrame(
            {
                'GCCC': [259, 0, 0],
                'AGTA': [0, 350, 0],
                'GGGC': [0, 0, 259],
                'TACT': [0, 0, 0]
            }, index=['GGGC', 'TACT', 'GCCC'])
        inputs.update({'fidelity_data': invalid_fidelity_data})
        with self.assertRaises(ValueError) as context:
            assign_fusion_sites(**inputs)
        self.assertEqual("fidelity_data.shape[0] != fidelity_data.shape[1]", str(context.exception))


class TestSelectJunctionByCodonUsage(TestCase):

    def setUp(self) -> None:
        self.junction_dna_map_sliding_window = [
            {'junction_dna': 'TGGTGGAAAAAA', 'i': 3, 'fusion_site': 'TGGA'},  # 1*1*0.39*0.39
            {'junction_dna': 'TGGTGGAAAAAG', 'i': 3, 'fusion_site': 'TGGA'},  # 1*1*0.39*0.61
            {'junction_dna': 'TGGTGGAAGAAA', 'i': 3, 'fusion_site': 'TGGA'},  # 1*1*0.61*0.39
            {'junction_dna': 'TGGTGGAAGAAG', 'i': 3, 'fusion_site': 'TGGA'},  # 1*1*0.61*0.61

        ]

    def test_select_junction_by_codon_usage_normal_cases(self):

        expected_sel_juction = {'junction_dna': 'TGGTGGAAGAAG', 'i': 3, 'fusion_site': 'TGGA'}
        sel_juction = select_junction_by_codon_usage(junctions=self.junction_dna_map_sliding_window,
                                                     codon_usage_table_path=codon_usage_table_path)

        self.assertEqual(expected_sel_juction, sel_juction)

    def test_select_junction_by_codon_usage_edge_cases(self):
        # edge cases: empty junction list

        edge_cases = [
            ('junctions', []),
        ]

        for param, value in edge_cases:

            inputs = dict(junctions=self.junction_dna_map_sliding_window,
                          codon_usage_table_path=codon_usage_table_path)
            inputs.update({param: value})
            expected_exceptions = {
                'junctions': f'No junctions are provided!'
            }

            if param in expected_exceptions.keys():
                with self.assertRaises(ValueError) as context:
                    select_junction_by_codon_usage(junctions=inputs['junctions'],
                                                   codon_usage_table_path=codon_usage_table_path)
                self.assertEqual(expected_exceptions[param], str(context.exception))


class TestConcatSelFusionSitesToFragments(TestCase):
    # example amino acid sequence is a part of published Human IgG1: https://rest.uniprot.org/uniprotkb/P0DOX5.fasta

    def setUp(self) -> None:
        self.fusion_sites = [
            "AGGA",
            "AAGA"
        ]
        self.fragments = [
            "QVQLVQSGG",
            "GVVQPGRSLR",
            "LSCAASGFTFSRYTIHWVRQA"
        ]
        self.junction_dna_map_sliding_window = [
            [{'junction_dna': 'GGAGGAGGAGTA', 'i': 5, 'fusion_site': 'AGGA'},
             {'junction_dna': 'GGAGGAGGAGTG', 'i': 5, 'fusion_site': 'AGGA'},
             {'junction_dna': 'GGGGGAGGAGTA', 'i': 5, 'fusion_site': 'AGGA'},
             {'junction_dna': 'GGGGGAGGAGTG', 'i': 5, 'fusion_site': 'AGGA'}],
            [{'junction_dna': 'TTAAGACTTAGA', 'i': 2, 'fusion_site': 'AAGA'},
             {'junction_dna': 'TTAAGACTTAGG', 'i': 2, 'fusion_site': 'AAGA'},
             {'junction_dna': 'TTAAGACTGAGA', 'i': 2, 'fusion_site': 'AAGA'},
             {'junction_dna': 'TTAAGACTGAGG', 'i': 2, 'fusion_site': 'AAGA'}]
        ]

        self.host = 'c_griseus'

    def test_concat_sel_fusion_sites_to_fragments_normal_cases(self):
        expected_fragment_with_fusion_sites = {
            "QVQLVQSGG": {
                "order": 0,
                "c_term_dna": "GG<A|GGA>",
                "middle_aa": "QVQLVQSG"
            },
            "GVVQPGRSLR": {
                "order": 1,
                "n_term_dna": "<A|GGA>",
                "middle_aa": "VVQPGRS",
                "c_term_dna": "TT<AAGA|>"
            },
            "LSCAASGFTFSRYTIHWVRQA": {
                "order": 2,
                "n_term_dna": "<AAGA|>",
                "middle_aa": "LSCAASGFTFSRYTIHWVRQA"
            }
        }

        fragment_with_fusion_sites = concat_sel_fusion_sites_to_fragments(
            fragments=self.fragments,
            fusion_sites=self.fusion_sites,
            sel_junction_dna_map_fusion_sites=self.junction_dna_map_sliding_window,
            codon_usage_table_path=codon_usage_table_path
        )

        self.assertEqual(expected_fragment_with_fusion_sites, fragment_with_fusion_sites)
