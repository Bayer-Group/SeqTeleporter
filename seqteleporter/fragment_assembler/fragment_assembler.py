from collections import defaultdict
import itertools
import random
import re
from Bio.Seq import Seq
from Bio.Restriction import *
from typing import List, Union

from seqteleporter.utils.utils import find_point_mutations, find_first_stop_codon_idx
from seqteleporter.config import ENZYME_INFO


def group_mutant_dna_frags_by_frag_span(mutant_dna_fragments: list) -> list:
    # Create a defaultdict with lists as default values
    grouped = defaultdict(list)

    # Group the dictionaries by the 'key' value
    for frag in mutant_dna_fragments:
        grouped[re.sub("_.+$", "", frag['name'])].append(frag)
    # Extract the grouped dictionaries into a list of lists
    result = list(grouped.values())

    return result


def generate_all_possible_variants_from_modules(module_names: list) -> list[list]:
    # Create a defaultdict with lists as default values
    grouped = defaultdict(list)
    # Group the dictionaries by the 'key' value
    for name in module_names:
        grouped[re.sub("_.+$", "", name)].append(name)
    # Extract the grouped dictionaries into a list of lists
    result = list(grouped.values())
    variants = []
    for frags in itertools.product(*result):
        sorted_frags = sorted(list(frags), key=lambda x: re.sub('-.+', '', x))
        sorted_frags = [re.sub('_', ',', re.sub('^.+?_', '', i)) for i in sorted_frags if not re.search('wild_type', i)]
        variant_name = ','.join(sorted_frags)
        muts_list = variant_name.split(',')
        if muts_list != ['']:
            variants.append(muts_list)
    return variants


def random_fragment_combinations(mutant_dna_fragments: list) -> list:
    grouped_dna_frags = group_mutant_dna_frags_by_frag_span(mutant_dna_fragments=mutant_dna_fragments)
    random_combi = [random.choice(lst) for lst in grouped_dna_frags]
    return random_combi


class AssembledFragment:
    def __init__(self, name, dna, aa, coding_start, coding_end, wt_aa, n_term_aa, c_term_aa):
        self.name = name
        self.dna = dna
        self.aa = aa
        self.n_term_aa = n_term_aa
        self.c_term_aa = c_term_aa
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.wt_aa = wt_aa
        self.translated_dna_assembly = str(Seq(self.dna)[self.coding_start:self.coding_end].translate())

    def validate_seq_length(self) -> bool:
        if len(self.wt_aa) == len(self.aa):
            return True
        print(f'validate_seq_length() failed: len(self.wt_aa)={len(self.wt_aa)}; len(self.aa)={len(self.aa)}')
        return False

    def validate_name_with_aa(self) -> bool:
        # the numbering of mutations does not consider the AAs encoded by the 5'DNA and 3'DNA, thus the offset.
        n_term_offset = len(self.n_term_aa)
        c_term_offset = len(self.c_term_aa)
        if c_term_offset == 0:
            point_mutations = find_point_mutations(wt_seq=self.wt_aa[n_term_offset:],
                                                   mut_seq=self.aa[n_term_offset:])
        else:
            point_mutations = find_point_mutations(wt_seq=self.wt_aa[n_term_offset:-c_term_offset],
                                                   mut_seq=self.aa[n_term_offset:-c_term_offset])
        mutations_in_name = re.findall(r"[A-Z]\d+[A-Z]", self.name)
        if set(point_mutations) == set(mutations_in_name):
            return True
        print(f'validate_name_with_aa() failed: '
              f'{set(point_mutations)-set(mutations_in_name)}; '
              f'{set(mutations_in_name)-set(point_mutations)}')
        print(self.name, self.wt_aa, self.aa)
        return False

    def validate_translated_dna_with_aa(self):
        difference = find_point_mutations(wt_seq=self.aa, mut_seq=self.translated_dna_assembly)
        if len(difference) == 0:
            return True
        print(f'validate_translated_dna_with_aa() failed: {difference}')
        return False

    def validate_assembled_fragment(self) -> bool:
        validated = \
            self.validate_seq_length() & \
            self.validate_name_with_aa() & \
            self.validate_translated_dna_with_aa()
        return validated


def assemble_fragments(a_set_of_dna_fragments: List[dict], enzyme: str, wt_seq: str,
                       five_prime_dna: str, three_prime_dna: str) -> AssembledFragment:

    # Find the AA encoded by 5'DNA and 3'DNA
    biopy_enzyme = globals()[enzyme]
    # because Biopython Restriction pyckage does not catalyse the reaction when there is less than 5 extra bases
    # adjacent to cutting site (this is to reflect in-vitro reaction), an "AAAAA" stuffer is temporarily added to allow
    # cutting simulation.
    stuffer = 'A' * 5
    five_prime_dna_digested_frags = biopy_enzyme.catalyse(Seq(five_prime_dna)+stuffer)
    if len(five_prime_dna_digested_frags) == 2:
        five_prime_dna_digested = re.sub(stuffer, '', str(five_prime_dna_digested_frags[1]))
    else:
        # enzyme site is not found or found too many times in 5'_dna
        raise ValueError(f'{enzyme} is not found or found too many times in DNA_5_PRIME({five_prime_dna})!')

    three_prime_dna_digested_frags = biopy_enzyme.catalyse(Seq(three_prime_dna))

    if len(three_prime_dna_digested_frags) == 2:
        three_prime_dna_digested = three_prime_dna_digested_frags[0]
    # Tackle edge case when the cut happens at the very edge
    elif len(three_prime_dna_digested_frags) == 1:
        # because Biopython Restriction pyckage does not catalyse the reaction when there is less than 5 extra bases
        # adjacent to cutting site (this is to reflect in-vitro reaction), an "AAAAA" stuffer is temporarily added to
        # allow cutting simulation.
        digested = biopy_enzyme.catalyse(stuffer + Seq(three_prime_dna))
        if digested[0] == stuffer:
            three_prime_dna_digested = Seq("")
        # possibility 1: The digestion happened at the very end of the sequence: ex: "TAATAGAGACCTTTAA". In this case,
        # add stuffer AAAAA to and then the digestion produces two fragments, with the first fragment == stuffer.

        # possibility 2: The digestion did not happen because the enzyme recognition site does not exist. In this case,
        # add stuffer AAAAA to and then the digested[0] != stuffer
        else:
            raise ValueError(f'DNA_3_PRIME({three_prime_dna}) is not digested by {enzyme}!')
    else:
        raise ValueError(f'DNA_3_PRIME({three_prime_dna}) is digested too many times by {enzyme}!')

    find_atg = re.search('ATG', str(five_prime_dna_digested))

    if find_atg:  # if find ATG in digested 5'DNA, that's the coding start.
        coding_start = find_atg.start()
    else:  # if no ATG in digested 5'DNA, coding starts after 5'DNA
        coding_start = len(five_prime_dna_digested)
    n_term_aa = str(Seq(five_prime_dna_digested[coding_start:]).translate())
    stop_codon_range = find_first_stop_codon_idx(three_prime_dna_digested)
    if stop_codon_range:
        c_term_aa = str(three_prime_dna_digested[:stop_codon_range[0]].translate())
        number_of_c_nt_after_stop_codon = len(three_prime_dna_digested) - stop_codon_range[0]
    else:
        c_term_aa = str(three_prime_dna_digested.translate())
        number_of_c_nt_after_stop_codon = 0

        # Validate that each fragment is cut exactly twice by the specified enzyme
    digested_frags = [biopy_enzyme.catalyse(Seq(frag['concatenated_dna'])) for frag in a_set_of_dna_fragments]
    for digested_frag in digested_frags:
        if len(digested_frag) != 3:
            raise ValueError(
                f'ERROR! Fragment is NOT cut exactly twice by the specified enzyme. Here are the digested fragment(s): '
                f'{digested_frag}')

    # Concat all digested fragments and store it as and object of class AssembledFragment
    dna = Seq('').join([frag[1] for frag in digested_frags])
    _name: list = []
    for frag in a_set_of_dna_fragments:
        _name = _name + re.findall(r"[A-Z]\d+[A-Z]", frag['name'])
    if len(_name) == 0:
        name = 'wild_type'
    else:
        name = '_'.join(_name)

    assembled_frag: AssembledFragment = AssembledFragment(
        name=name,
        dna=dna,
        aa="".join([n_term_aa]+[f['fragment'] for f in a_set_of_dna_fragments]+[c_term_aa]),
        n_term_aa=n_term_aa,
        c_term_aa=c_term_aa,
        coding_start=coding_start,
        coding_end=len(dna)-number_of_c_nt_after_stop_codon,
        wt_aa=''.join([n_term_aa, wt_seq, c_term_aa])
    )
    return assembled_frag
