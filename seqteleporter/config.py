from typing import Dict, Union
from seqteleporter.utils.utils import expand_python_codon_tables


# EXPAND HOSTS FOR EXTERNAL PACKAGE
addtional_hosts = [('p_pastoris', 4922), ('c_griseus', 10029)]
for host in addtional_hosts:
    expand_python_codon_tables(*host)

PLATE_FORMATS = {
    24: {'rows': ('A', 'D'), 'columns': (1, 6)},
    80: {'rows': ('A', 'H'), 'columns': (1, 10)},
    96: {'rows': ('A', 'H'), 'columns': (1, 12)},
    384: {'rows': ('A', 'P'), 'columns': (1, 24)},
}

ENZYME_INFO: Dict[str, Union[Dict[str, Union[str, int]], Dict[str, Union[str, int]]]] = {
    'BsaI': {'fusion_site_length': 4,
             'recognition_site': 'GGTCTCN'},
    'SapI': {'fusion_site_length': 3,
             'recognition_site': 'GCTCTTCN'},
    'BpiI': {'fusion_site_length': 4,
             'recognition_site': 'GAAGACNN'}
}

PARTITION_SEARCH_MODES = {
    'cost_scan': {'cost_scan': True, 'pre_distribute_mutations': False, 'one_dist': False},
    'dist_mut': {'cost_scan': False, 'pre_distribute_mutations': True, 'one_dist': False},
    'dist_mut_cost_scan': {'cost_scan': True, 'pre_distribute_mutations': True, 'one_dist': False},
    'dist_mut_1': {'cost_scan': False, 'pre_distribute_mutations': True, 'one_dist': True},
    'dist_mut_1_cost_scan': {'cost_scan': True, 'pre_distribute_mutations': True, 'one_dist': True},
    'exhaustive': {'cost_scan': False, 'pre_distribute_mutations': False, 'one_dist': False},
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


class DnaProviderSpec:
    def __init__(self, provider_name: str, product_name: str, min_len: int, max_len: int):
        self.provider_name = provider_name
        self.product_name = product_name
        self.min_len = min_len
        self.max_len = max_len


EBLOCK = DnaProviderSpec(
    provider_name='IDT',
    product_name='eBlock',
    min_len=300,
    max_len=1500,
)

GBLOCK = DnaProviderSpec(
    provider_name='IDT',
    product_name='gBlock',
    min_len=125,
    max_len=3000,
)

GENEART = DnaProviderSpec(
    provider_name='GeneArt',
    product_name='GeneArt Sub-cloned',
    min_len=200,
    max_len=12000,
)

