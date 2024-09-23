"""Load Input Parameters"""

import ast
import os
from os.path import join, dirname, abspath
from typing import Dict, Union

from dotenv import load_dotenv

print('--------------------------------------------------------------------------------------------------')
print(f'Load input_draft parameters.')
print('--------------------------------------------------------------------------------------------------')

CONFIG_DIR = dirname(abspath(__file__))
PROJECT_ROOT = dirname(dirname(CONFIG_DIR))
INPUT_PATH = join(CONFIG_DIR, 'partition_search_opt_input50.txt')
load_dotenv(INPUT_PATH)

ALLOWED_CUT_POSITIONS = ast.literal_eval(os.environ.get("ALLOWED_CUT_POSITIONS"))
print(f'Load ALLOWED_CUT_POSITIONS={ALLOWED_CUT_POSITIONS}')
DNA_5_PRIME = os.environ.get("DNA_5_PRIME")
print(f'Load DNA_5_PRIME={DNA_5_PRIME}')
DNA_3_PRIME = os.environ.get("DNA_3_PRIME")
print(f'Load DNA_3_PRIME={DNA_3_PRIME}')
FIX_DNA_SEQUENCE = os.environ.get("FIX_DNA_SEQUENCE")
print(f'Load FIX_DNA_SEQUENCE={FIX_DNA_SEQUENCE}')
FUSION_SITES_USED_BY_BACKBONE = ast.literal_eval(os.environ.get("FUSION_SITES_USED_BY_BACKBONE"))
HOST = os.environ.get("HOST")
print(f'Load HOST={HOST}')
FIDELITY_DATA = os.environ.get("FIDELITY_DATA")
FIDELITY_DATA_PATH = join(PROJECT_ROOT, 'data', 'neb_fidelity_data', FIDELITY_DATA)
print(f'Load FIDELITY_DATA_PATH={FIDELITY_DATA_PATH}')
SEQUENCE = os.environ.get("SEQUENCE")
print(f'Load SEQUENCE={SEQUENCE}')
MUTATIONS = ast.literal_eval(os.environ.get("MUTATIONS"))
print(f'Load MUTATIONS={MUTATIONS}')
LINKED_MUTATIONS = ast.literal_eval(os.environ.get("LINKED_MUTATIONS"))
print(f'Load LINKED_MUTATIONS={LINKED_MUTATIONS}')
CUT_NUMBER_RANGE = ast.literal_eval(os.environ.get("CUT_NUMBER_RANGE"))
print(f'Load CUT_NUMBER_RANGE={CUT_NUMBER_RANGE}')
MIN_FRAGMENT_LENGTH = int(os.environ.get("MIN_FRAGMENT_LENGTH"))
print(f'Load MIN_FRAGMENT_LENGTH={MIN_FRAGMENT_LENGTH}')
MAX_COST = float(os.environ.get("MAX_COST"))
print(f'Load MAX_COST={MAX_COST}')
MAX_LENGTH_UNEVENNESS = float(os.environ.get("MAX_LENGTH_UNEVENNESS"))
print(f'Load MAX_LENGTH_UNEVENNESS={MAX_LENGTH_UNEVENNESS}')
MIN_LIGATION_FIDELITY = float(os.environ.get("MIN_LIGATION_FIDELITY"))
print(f'Load MIN_LIGATION_FIDELITY={MIN_LIGATION_FIDELITY}')
SATISFACTION_LIGATION_FIDELITY = float(os.environ.get("SATISFACTION_LIGATION_FIDELITY"))
print(f'Load SATISFACTION_LIGATION_FIDELITY={SATISFACTION_LIGATION_FIDELITY}')

ENZYME_INFO: Dict[str, Union[Dict[str, Union[str, int]], Dict[str, Union[str, int]]]] = {
    'BsaI': {'fusion_site_length': 4,
             'recognition_site': 'GGTCTCN'},
    'SapI': {'fusion_site_length': 3,
             'recognition_site': 'GCTCTTCN'}
}

CODON_TABLE = {
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
