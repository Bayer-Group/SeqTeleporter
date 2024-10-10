from unittest import TestCase
import re
import warnings
from os.path import join, dirname, abspath
from os import listdir
from python_codon_tables.python_codon_tables import _tables_dir as codon_usage_tbl_dir

from seqteleporter.post_partition_processor.post_partition_processor import valid_end_dna, append_end_dna, \
    append_enzyme_sites_and_stuffer, make_mutant_aa_fragments, make_mutant_dna_fragments_from_mutant_aa_fragments
from seqteleporter.utils.utils import find_best_codon_by_usage

FIDELITY_DATA_PATH = join(
    dirname(dirname(dirname(abspath(__file__)))),
    'seqteleporter', 'data', 'neb_fidelity_data', 'FileS01_T4_01h_25C.xlsx'
)

HOST = 'e_coli'
for f in listdir(codon_usage_tbl_dir):
    if re.search(HOST, f):
        codon_usage_table_path = join(codon_usage_tbl_dir, f)
        break

# example amino acid sequence is a part of published Human IgG1: https://rest.uniprot.org/uniprotkb/P0DOX5.fasta
FRAGMENT_WITH_FUSION_SITES = {
    "QVQLVQSGG": {
        "order": 0,
        "n_term_dna": "(TTTGGTCTCCATGAAAGGTTCGCACCTG)",
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
        "middle_aa": "LSCAASGFTFSRYTIHWVRQA",
        "c_term_dna": "(CACACATCCCTGTCTGAGACCTTT)"
    }
}
MUTANT_AA_FRAGMENTS = [{'name': '1-9_wild_type', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QVQLVQSGG',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                        'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                        'middle_aa': 'QVQLVQSG'},
                       {'name': '1-9_V2P_S7A', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QPQLVQAGG',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                        'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                        'middle_aa': 'QPQLVQAG'},
                       {'name': '1-9_V5Q', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QVQLQQSGG',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                        'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                        'middle_aa': 'QVQLQQSG'},
                       {'name': '1-9_V2P_V5Q_S7A', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QPQLQQAGG',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                        'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                        'middle_aa': 'QPQLQQAG'},
                       {'name': '1-9_V5S', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QVQLSQSGG',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                        'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                        'middle_aa': 'QVQLSQSG'},
                       {'name': '1-9_V2P_V5S_S7A', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QPQLSQAGG',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                        'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                        'middle_aa': 'QPQLSQAG'},
                       {'name': '10-19_wild_type', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQPGRSLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'VVQPGRS'},
                       {'name': '10-19_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQPGRKLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'VVQPGRK'},
                       {'name': '10-19_P14T', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQTGRSLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'VVQTGRS'},
                       {'name': '10-19_P14T_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQTGRKLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'VVQTGRK'},
                       {'name': '10-19_P14M', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQMGRSLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'VVQMGRS'},
                       {'name': '10-19_P14M_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQMGRKLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'VVQMGRK'},
                       {'name': '10-19_V11I', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQPGRSLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'IVQPGRS'},
                       {'name': '10-19_V11I_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQPGRKLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'IVQPGRK'},
                       {'name': '10-19_V11I_P14T', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQTGRSLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'IVQTGRS'},
                       {'name': '10-19_V11I_P14T_S17K', 'fragment_start': 9, 'fragment_end': 19,
                        'fragment': 'GIVQTGRKLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'IVQTGRK'},
                       {'name': '10-19_V11I_P14M', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQMGRSLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'IVQMGRS'},
                       {'name': '10-19_V11I_P14M_S17K', 'fragment_start': 9, 'fragment_end': 19,
                        'fragment': 'GIVQMGRKLR',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                        'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)',
                        'middle_aa': 'IVQMGRK'},
                       {'name': '20-40_wild_type', 'fragment_start': 19, 'fragment_end': 40,
                        'fragment': 'LSCAASGFTFSRYTIHWVRQA',
                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTTTGGTCTCT)<AAGA|>',
                        'c_term_dna': '(CACACATCCCTGTCTGAGACCTTT)(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGT)',
                        'middle_aa': 'LSCAASGFTFSRYTIHWVRQA'}]
FRAGMENTS_N_AND_C_TERM_DNA = {'QVQLVQSGG': {'order': 0,
                                            'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                                            'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)',
                                            'middle_aa': 'QVQLVQSG'},
                              'GVVQPGRSLR': {'order': 1,
                                             'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>',
                                             'middle_aa': 'VVQPGRS',
                                             'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)'},
                              'LSCAASGFTFSRYTIHWVRQA': {'order': 2,
                                                        'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTTTGGTCTCT)<AAGA|>',
                                                        'middle_aa': 'LSCAASGFTFSRYTIHWVRQA',
                                                        'c_term_dna': '(CACACATCCCTGTCTGAGACCTTT)(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGT)'}}
MUTANT_DNA_FRAGMENTS = [
    {'name': '1-9_wild_type', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QVQLVQSGG', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)', 'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)', 'middle_aa': 'QVQLVQSG', 'middle_dna': 'CAGGTGCAGCTGGTGCAGAGCGGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTTTTGGTCTCCATGAAAGGTTCGCACCTGCAGGTGCAGCTGGTGCAGAGCGGCGGAGGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT'},
    {'name': '1-9_V2P_S7A', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QPQLVQAGG', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)', 'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)', 'middle_aa': 'QPQLVQAG', 'middle_dna': 'CAGCCGCAGCTGGTGCAGGCGGGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTTTTGGTCTCCATGAAAGGTTCGCACCTGCAGCCGCAGCTGGTGCAGGCGGGCGGAGGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT'},
    {'name': '1-9_V5Q', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QVQLQQSGG', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)', 'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)', 'middle_aa': 'QVQLQQSG', 'middle_dna': 'CAGGTGCAGCTGCAGCAGAGCGGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTTTTGGTCTCCATGAAAGGTTCGCACCTGCAGGTGCAGCTGCAGCAGAGCGGCGGAGGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT'},
    {'name': '1-9_V2P_V5Q_S7A', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QPQLQQAGG', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)', 'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)', 'middle_aa': 'QPQLQQAG', 'middle_dna': 'CAGCCGCAGCTGCAGCAGGCGGGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTTTTGGTCTCCATGAAAGGTTCGCACCTGCAGCCGCAGCTGCAGCAGGCGGGCGGAGGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT'},
    {'name': '1-9_V5S', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QVQLSQSGG', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)', 'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)', 'middle_aa': 'QVQLSQSG', 'middle_dna': 'CAGGTGCAGCTGAGCCAGAGCGGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTTTTGGTCTCCATGAAAGGTTCGCACCTGCAGGTGCAGCTGAGCCAGAGCGGCGGAGGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT'},
    {'name': '1-9_V2P_V5S_S7A', 'fragment_start': 0, 'fragment_end': 9, 'fragment': 'QPQLSQAGG', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)', 'c_term_dna': 'GG<A|GGA>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT)', 'middle_aa': 'QPQLSQAG', 'middle_dna': 'CAGCCGCAGCTGAGCCAGGCGGGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTTTTGGTCTCCATGAAAGGTTCGCACCTGCAGCCGCAGCTGAGCCAGGCGGGCGGAGGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCT'},
    {'name': '10-19_wild_type', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQPGRSLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'VVQPGRS', 'middle_dna': 'GTGGTGCAGCCGGGCCGCAGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAGTGGTGCAGCCGGGCCGCAGCTTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQPGRKLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'VVQPGRK', 'middle_dna': 'GTGGTGCAGCCGGGCCGCAAA', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAGTGGTGCAGCCGGGCCGCAAATTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_P14T', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQTGRSLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'VVQTGRS', 'middle_dna': 'GTGGTGCAGACCGGCCGCAGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAGTGGTGCAGACCGGCCGCAGCTTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_P14T_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQTGRKLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'VVQTGRK', 'middle_dna': 'GTGGTGCAGACCGGCCGCAAA', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAGTGGTGCAGACCGGCCGCAAATTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_P14M', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQMGRSLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'VVQMGRS', 'middle_dna': 'GTGGTGCAGATGGGCCGCAGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAGTGGTGCAGATGGGCCGCAGCTTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_P14M_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GVVQMGRKLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'VVQMGRK', 'middle_dna': 'GTGGTGCAGATGGGCCGCAAA', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAGTGGTGCAGATGGGCCGCAAATTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_V11I', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQPGRSLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'IVQPGRS', 'middle_dna': 'ATTGTGCAGCCGGGCCGCAGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAATTGTGCAGCCGGGCCGCAGCTTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_V11I_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQPGRKLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'IVQPGRK', 'middle_dna': 'ATTGTGCAGCCGGGCCGCAAA', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAATTGTGCAGCCGGGCCGCAAATTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_V11I_P14T', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQTGRSLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'IVQTGRS', 'middle_dna': 'ATTGTGCAGACCGGCCGCAGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAATTGTGCAGACCGGCCGCAGCTTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_V11I_P14T_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQTGRKLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'IVQTGRK', 'middle_dna': 'ATTGTGCAGACCGGCCGCAAA', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAATTGTGCAGACCGGCCGCAAATTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_V11I_P14M', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQMGRSLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'IVQMGRS', 'middle_dna': 'ATTGTGCAGATGGGCCGCAGC', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAATTGTGCAGATGGGCCGCAGCTTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '10-19_V11I_P14M_S17K', 'fragment_start': 9, 'fragment_end': 19, 'fragment': 'GIVQMGRKLR', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCT)<A|GGA>', 'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA)', 'middle_aa': 'IVQMGRK', 'middle_dna': 'ATTGTGCAGATGGGCCGCAAA', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTTTTGGTCTCTAGGAATTGTGCAGATGGGCCGCAAATTAAGATGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTAGACGTAGTCGGACCTGCATAGGTA'},
    {'name': '20-40_wild_type', 'fragment_start': 19, 'fragment_end': 40, 'fragment': 'LSCAASGFTFSRYTIHWVRQA', 'n_term_dna': '(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTTTGGTCTCT)<AAGA|>', 'c_term_dna': '(CACACATCCCTGTCTGAGACCTTT)(ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGT)', 'middle_aa': 'LSCAASGFTFSRYTIHWVRQA', 'middle_dna': 'CTGAGCTGCGCGGCGAGCGGCTTTACCTTTAGCCGCTATACCATTCATTGGGTGCGCCAGGCG', 'concatenated_dna': 'ATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGTTTGGTCTCTAAGACTGAGCTGCGCGGCGAGCGGCTTTACCTTTAGCCGCTATACCATTCATTGGGTGCGCCAGGCGCACACATCCCTGTCTGAGACCTTTATCGCGGGTTTTTTCCCTGAGATCCCATCAGATAAATTGGAGTATCCCGTGCACGGACTGCTCCCGATAGGATTTCACTCTCCGGAATTTGAGGGACAGT'}
]

class TestValidEndDna(TestCase):

    def test_valid_end_dna_2_warnings(self):
        five_prime_dna = "TTTGTCTCCATGAAAGGTTCGCACCTG"
        three_prime_dna = "CACACATCCCTGTCTGAGACCTTT"
        enzyme = "BsaI"
        fusion_sites_used_by_backbone = ('CCAT', 'TAAT')
        with warnings.catch_warnings(record=True) as w:
            valid_end_dna(five_prime_dna, three_prime_dna, enzyme, fusion_sites_used_by_backbone)
            assert str(w[0].message) == f"Warning! {enzyme} site is not found in 5'-sequence."
            assert str(w[1].message) == (f"Warning! The overhang produced by {enzyme} digestion of 3'-sequence, "
                                         f"i.e. TCCC), is not found in the specified fusion sites used by backbone, "
                                         f"i.e., {fusion_sites_used_by_backbone}.")


class TestAppendEndDna(TestCase):

    def test_append_end_dna(self):
        five_prime_dna = "TTTGGTCTCCATGAAAGGTTCGCACCTG"
        three_prime_dna = "CACACATCCCTGTCTGAGACCTTT"
        fragment_with_fusion_sites = {
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
        expected_output = FRAGMENT_WITH_FUSION_SITES
        output = append_end_dna(fragment_with_fusion_sites, five_prime_dna, three_prime_dna)
        self.assertEqual(expected_output, output)


class TestAppendEnzymeSitesAndStuffer(TestCase):

    def test_append_enzyme_sites_and_stuffer(self):
        enzyme = "BsaI"
        min_dna_frag_length = 300
        expected_output = {'QVQLVQSGG': {'order': 0,
                                         'n_term_dna': '(AGGTAACCCAGTTCTTAGCACACATCCGTTTTCTCTATGACCACGCTCGATGTCGATCGCCTCAATTAGCGGACTTGTGTGCGTTAATGTGCTCCGTTGGGTTGCCCCCAAGAAGT)(TTTGGTCTCCATGAAAGGTTCGCACCTG)',
                                         'c_term_dna': 'GG<A|GGA>(TGAGACCTTTCGCCAAGAATTCATCGTAAGTGACCTGCGCTGTGGGCGACTATAGTGGTTTCTTCTGACGCAGCCTACCTGCGGTTTGAATTGCTGGTCACATCGCTCTTCGACTTGCGGATGAAA)',
                                         'middle_aa': 'QVQLVQSG'}, 'GVVQPGRSLR': {'order': 1,
                                                                                  'n_term_dna': '(AGGTAACCCAGTTCTTAGCACACATCCGTTTTCTCTATGACCACGCTCGATGTCGATCGCCTCAATTAGCGGACTTGTGTGCGTTAATGTGCTCCGTTGGGTTGCCCCCAAGAAGTCGCCAAGATTTGGTCTCT)<A|GGA>',
                                                                                  'middle_aa': 'VVQPGRS',
                                                                                  'c_term_dna': 'TT<AAGA|>(TGAGACCTTTATTCATCGTAAGTGACCTGCGCTGTGGGCGACTATAGTGGTTTCTTCTGACGCAGCCTACCTGCGGTTTGAATTGCTGGTCACATCGCTCTTCGACTTGCGGATGAAAACGCTTGCAGGCGAGCT)'},
                           'LSCAASGFTFSRYTIHWVRQA': {'order': 2,
                                                     'n_term_dna': '(AGGTAACCCAGTTCTTAGCACACATCCGTTTTCTCTATGACCACGCTCGATGTCGATCGCCTCAATTAGCGGACTTGTGTGCGTTAATGTGCTCCGTTGTTTGGTCTCT)<AAGA|>',
                                                     'middle_aa': 'LSCAASGFTFSRYTIHWVRQA',
                                                     'c_term_dna': '(CACACATCCCTGTCTGAGACCTTT)(GGTTGCCCCCAAGAAGTCGCCAAGAATTCATCGTAAGTGACCTGCGCTGTGGGCGACTATAGTGGTTTCTTCTGACGCAGCCTACCTGCGGTTTGAATTG)'}}
        output = append_enzyme_sites_and_stuffer(FRAGMENT_WITH_FUSION_SITES, enzyme, min_dna_frag_length)
        self.assertEqual(expected_output, output)
        patterns = r"[()<>|]"
        for frag in output.values():
            current_nterm_frag_dna_length = len(re.sub(patterns, "", frag.get("n_term_dna", "")))
            current_cterm_frag_dna_length = len(re.sub(patterns, "", frag.get("c_term_dna", "")))
            current_middle_frag_dna_length = len(frag.get("middle_aa")) * 3
            current_frag_dna_length = sum([current_nterm_frag_dna_length,
                                           current_cterm_frag_dna_length,
                                           current_middle_frag_dna_length])
            # print(f'append_enzyme_sites_and_stuffer: {current_frag_dna_length}')
            self.assertEqual(min_dna_frag_length, current_frag_dna_length)


class TestMakeFragmentsWithMutations(TestCase):

    def test_make_fragments_with_mutations(self):
        mutations_0idx = [{'position': 4, 'aa': ['Q', 'S']},
                          # {'position': 8, 'aa': ['A']},  # should return error Can not generate mutation that fits the constraint of fusion site!
                          {'position': 10, 'aa': ['I']},
                          {'position': 13, 'aa': ['T', 'M']},
                          {'position': 16, 'aa': ['K']}]
        linked_mutations_0idx = [(('V', 1, 'P'), ('S', 6, 'A'))]
        expected_output = MUTANT_AA_FRAGMENTS
        output = make_mutant_aa_fragments(FRAGMENTS_N_AND_C_TERM_DNA, mutations_0idx, linked_mutations_0idx,
                                          codon_usage_table_path, positions_include_wt_aa_0idx=[1, 4, 6, 10, 13, 16])

        self.assertEqual(sorted(expected_output, key=lambda d: d['name']), sorted(output, key=lambda d: d['name']))


class TestMakeDnaFragmentsFromAaFragments(TestCase):
    def test_make_dna_fragments_from_aa_fragments(self):
        aas = "".join([k for k in FRAGMENTS_N_AND_C_TERM_DNA.keys()])
        fix_base = {0: None, 1: None, 2: None}
        specify_dna_seq = "".join([find_best_codon_by_usage(codon_usage_table_path, aa, fix_base) for aa in aas])
        expected_output = MUTANT_DNA_FRAGMENTS
        output = make_mutant_dna_fragments_from_mutant_aa_fragments(MUTANT_AA_FRAGMENTS,
                                                                    codon_usage_table_path,
                                                                    specify_dna_seq, "BsaI")

        self.assertEqual(expected_output, output)
