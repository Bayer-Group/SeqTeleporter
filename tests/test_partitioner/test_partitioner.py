import shutil
from os import mkdir
from os.path import dirname, abspath, join, exists
from unittest import TestCase

from proseqteleporter.config import CODON_TABLE, PARTITION_SEARCH_MODES
from proseqteleporter.partitioner.compute_best_partitions import compute_best_partitions


FIDELITY_DATA_PATH = join(
    dirname(dirname(dirname(abspath(__file__)))),
    'proseqteleporter','data','neb_fidelity_data','FileS01_T4_01h_25C.xlsx'
)


class TestComputeBestPartitions(TestCase):

    # example amino acid sequence is a part of published Human IgG1: https://rest.uniprot.org/uniprotkb/P0DOX5.fasta
    def setUp(self) -> None:

        self.output_dir = join(dirname(abspath(__file__)), 'output_secret')
        self.log_dir = join(self.output_dir, 'logs')
        self.result_dir = join(self.output_dir, 'results')
        self.inputs = dict(
            s="QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQA",
            mutations_0idx=[{'position': 4, 'aa': ['Q', 'S']},
                            {'position': 13, 'aa': ['T', 'M']},
                            {'position': 17, 'aa': ['T', 'M']},
                            {'position': 20, 'aa': ['T', 'M']},
                            {'position': 34, 'aa': ['G']}],
            linked_mutations_0idx=[(('Q', 5, 'K'), ('G', 8, 'H'))],
            cut_number_range=(4, 5),
            fidelity_data_path=FIDELITY_DATA_PATH,
            min_aa_length=6,
            fusion_sites_used_by_backbone=('CCAT', 'TGTC', 'TAAT'),
            max_cost=2000,
            max_unevenness=1,
            min_ligation_fidelity=0.9,
            satisfaction_fidelity=0.95,
            supress_output=False,
            search_method="BFS",
            codon_usage_tbl_dir=join(dirname(dirname(FIDELITY_DATA_PATH)), 'codon_usage'),
            host='c_griseus',
            sort_by_cost=False,
            enzyme='BsaI',
            allowed_cut_positions_1idx=[],
            partition_search_mode='exhaustive',
            select_top_n_partitions=3,
            cost_per_nt=0.06,
            provider_min_frag_len=300,
            provider_max_frag_len=1500,
            output_dir=join(dirname(abspath(__file__)), 'output_secret'),
            max_partition_number_checked=100000
        )
        # create output_secret dirs if not exist
        for dir_ in [self.output_dir, self.log_dir, self.result_dir]:
            if not exists(dir_):
                mkdir(dir_)

    def test_compute_best_partitions_normal_case(self):
        expected_output = {
            'sequence': 'QVQLVQSGGGVVQPGRSLRLSCAASGFTFSRYTIHWVRQA',
            'mutations': [{'position': 5, 'aa': ['Q', 'S']}, {'position': 14, 'aa': ['T', 'M']}, {'position': 18, 'aa': ['T', 'M']}, {'position': 21, 'aa': ['T', 'M']}, {'position': 35, 'aa': ['G']}],
            'diversity': 324,
            'best_partitions_by_cut_number': [
                {'number_of_cuts': 4,
                 'number_of_partitions_checked': 3876,
                 'sel_partitions': [
                     {'partition': (9, 16, 23, 32),
                      'ligation_fidelity': 0.99,
                      'fragment_length_unevenness': 0.286,
                      'cost': 324.0,
                      'fusion_sites': ('AGGA', 'AAGA', 'AGCA', 'ACAA'),
                      'fragments': ['QVQLVQSGG', 'GVVQPGR', 'SLRLSCA', 'ASGFTFSRY', 'TIHWVRQA'],
                      'expression': '\n[QVQLVQSG]GG<A|GGA>\n            <A|GGA>[VVQP]GG<AAGA|>\n                           <AAGA|>[SLRLSC]GC<A|GCA>\n                                            <A|GCA>[SGFTFSRY]<|ACAA>\n                                                             <|ACAA>TA[HWVRQA]',
                      'fragment_with_fusion_sites': {'QVQLVQSGG': {'order': 0, 'c_term_dna': 'GG<A|GGA>', 'middle_aa': 'QVQLVQSG'}, 'GVVQPGR': {'order': 1, 'n_term_dna': '<A|GGA>', 'middle_aa': 'VVQP', 'c_term_dna': 'GG<AAGA|>'}, 'SLRLSCA': {'order': 2, 'n_term_dna': '<AAGA|>', 'middle_aa': 'SLRLSC', 'c_term_dna': 'GC<A|GCA>'}, 'ASGFTFSRY': {'order': 3, 'n_term_dna': '<A|GCA>', 'middle_aa': 'SGFTFSRY', 'c_term_dna': '<|ACAA>'}, 'TIHWVRQA': {'order': 4, 'n_term_dna': '<|ACAA>TA', 'middle_aa': 'HWVRQA'}}},
                     {'partition': (9, 16, 24, 31),
                      'ligation_fidelity': 0.99,
                      'fragment_length_unevenness': 0.286,
                      'cost': 324.0,
                      'fusion_sites': ('AGGA', 'AAGA', 'AAGC', 'AGAT'),
                      'fragments': ['QVQLVQSGG', 'GVVQPGR', 'SLRLSCAA', 'SGFTFSR', 'YTIHWVRQA'],
                      'expression': '\n[QVQLVQSG]GG<A|GGA>\n            <A|GGA>[VVQP]GG<AAGA|>\n                           <AAGA|>[SLRLSCA]GC<A|AGC>\n                                             <A|AGC>[GFTFS]<AGA|T>\n                                                           <AGA|T>AC[TIHWVRQA]', 'fragment_with_fusion_sites': {'QVQLVQSGG': {'order': 0, 'c_term_dna': 'GG<A|GGA>', 'middle_aa': 'QVQLVQSG'}, 'GVVQPGR': {'order': 1, 'n_term_dna': '<A|GGA>', 'middle_aa': 'VVQP', 'c_term_dna': 'GG<AAGA|>'}, 'SLRLSCAA': {'order': 2, 'n_term_dna': '<AAGA|>', 'middle_aa': 'SLRLSCA', 'c_term_dna': 'GC<A|AGC>'}, 'SGFTFSR': {'order': 3, 'n_term_dna': '<A|AGC>', 'middle_aa': 'GFTFS', 'c_term_dna': '<AGA|T>'}, 'YTIHWVRQA': {'order': 4, 'n_term_dna': '<AGA|T>AC', 'middle_aa': 'TIHWVRQA'}}},
                     {'partition': (9, 16, 24, 32),
                      'ligation_fidelity': 0.99,
                      'fragment_length_unevenness': 0.286,
                      'cost': 324.0,
                      'fusion_sites': ('AGGA', 'AAGA', 'AAGC', 'ACAA'),
                      'fragments': ['QVQLVQSGG', 'GVVQPGR', 'SLRLSCAA', 'SGFTFSRY', 'TIHWVRQA'],
                      'expression': '\n[QVQLVQSG]GG<A|GGA>\n            <A|GGA>[VVQP]GG<AAGA|>\n                           <AAGA|>[SLRLSCA]GC<A|AGC>\n                                             <A|AGC>[GFTFSRY]<|ACAA>\n                                                             <|ACAA>TA[HWVRQA]', 'fragment_with_fusion_sites': {'QVQLVQSGG': {'order': 0, 'c_term_dna': 'GG<A|GGA>', 'middle_aa': 'QVQLVQSG'}, 'GVVQPGR': {'order': 1, 'n_term_dna': '<A|GGA>', 'middle_aa': 'VVQP', 'c_term_dna': 'GG<AAGA|>'}, 'SLRLSCAA': {'order': 2, 'n_term_dna': '<AAGA|>', 'middle_aa': 'SLRLSCA', 'c_term_dna': 'GC<A|AGC>'}, 'SGFTFSRY': {'order': 3, 'n_term_dna': '<A|AGC>', 'middle_aa': 'GFTFSRY', 'c_term_dna': '<|ACAA>'}, 'TIHWVRQA': {'order': 4, 'n_term_dna': '<|ACAA>TA', 'middle_aa': 'HWVRQA'}}}
                 ],
                 'num_of_checked_unique_partitions': 3876,
                 'hard_constraint_violations': {'Fragment(s) too short.': 3789}}
            ]
        }
        output, output_path_to_return = compute_best_partitions(**self.inputs)
        output.pop('total_elapsed_time')
        output['best_partitions_by_cut_number'][0].pop('elapsed_time')
        self.assertEqual(expected_output, output)
        shutil.rmtree(self.output_dir)

    def test_compute_best_partitions_invalid_input(self):

        invalid_inputs = (
            ('s', 'ABCDEF',
             'The provided input sequence is not a valid amino acid sequence!'),

            ('fusion_sites_used_by_backbone', ('CCAT', 'TGTM'),
             f"The fusion sites used by backbone are invalid! "
             f"Some of the fusion sites are not DNA sequence."
             ),

            ('fusion_sites_used_by_backbone', ('CCAT', 'ATGG', 'ATAG'),
             f"The fusion sites used by backbone are invalid! "
             f"Some of the fusion sites (or its reverse complement) used by backbone is duplicated."
             f"Please try again with valid fusion sites."
             ),
        )

        for param, val, errmsg in invalid_inputs:
            inputs = self.inputs.copy()
            inputs.update({param: val})
            with self.assertRaises(ValueError) as context:
                compute_best_partitions(**inputs)
            self.assertEqual(errmsg, str(context.exception))
