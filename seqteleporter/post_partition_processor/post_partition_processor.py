import warnings
import json
from os import path
from itertools import product
import re
import copy
from typing import Tuple, List, Union
import pandas as pd
from Bio.Restriction import *
from Bio.Seq import Seq, MutableSeq
from dnachisel import (DnaOptimizationProblem, AvoidPattern, AvoidHairpins, EnforceGCContent, CodonOptimize,
                       random_dna_sequence)

from seqteleporter.config import PLATE_FORMATS, ENZYME_INFO
from seqteleporter.utils.utils import (make_rev_compliment, find_best_codon_by_usage, generate_aa2codon_dict,
                                       remove_linked_mutations_from_all_mutations,
                                       multi_well_plate_position_generator, prepare_0idx_mutations)
from seqteleporter.utils.idt_tools import validate_seq_complexity_idt
from seqteleporter.utils.idt_config import IdtCredentials
from seqteleporter.utils.load_input_params import load_input_params
from seqteleporter.fragment_assembler.plate_mapper import load_module_sheet
from seqteleporter.fragment_assembler.fragment_assembler import random_fragment_combinations, assemble_fragments
from seqteleporter.partition_property_finder.cost_finder import find_exact_cost


# generated using generate_stuffer_sequence(seq_length=300)
STUFFER_REF = 'AGGTAACCCAGTTCTTAGCACACATCCGTTTTCTCTATGACCACGCTCGATGTCGATCGCCTCAATTAGCGGACTTGTGTGCGTTAATGTGCTCCGTTGGGT' \
              'TGCCCCCAAGAAGTCGCCAAGAATTCATCGTAAGTGACCTGCGCTGTGGGCGACTATAGTGGTTTCTTCTGACGCAGCCTACCTGCGGTTTGAATTGCTGGT' \
              'CACATCGCTCTTCGACTTGCGGATGAAAACGCTTGCAGGCGAGCTTTACACTCCTAATATAAGCGGGCGAACCTGGCACGAGATACACTCTCGTTACCCAGT' \
              'GTGGTTGTGCTGCTGTACAAGACCATACGCAACTTAGGATGCGGGGTTTTTCCAATGCCCTACAACGGGGTCGCTATATGGACTTTCCAAGGCGAGCCCCGT' \
              'TGATTTTATAAGTAGCGGCCTAAGTGAATGCGACGTAGCGGCGTAAGGAGGAACTATTTGTCTACGCACGCTCACCGTATACTCCAGATCCG'


def generate_stuffer_sequence(seq_length: int, enzyme: str) -> str:
    # Sequences that cause problems for gBlocks Gene Fragment synthesis are those having extremely low or high GC
    # content (less than 25% and greater than 75%), homopolymeric runs of 10 or more As and Ts or 6 or more Gs and
    # Cs. Other structural motifs such as repeats or hairpins may also influence our ability to construct gBlocks
    # Gene Fragments. The decision to accept or reject an order for a particular gBlocks Gene Fragment sequence is
    # complex and multifactorial. (source: IDT FAQs)

    # DEFINE THE OPTIMIZATION PROBLEM
    problem = DnaOptimizationProblem(
        sequence=random_dna_sequence(seq_length),
        constraints=[
            AvoidPattern("".join([enzyme, "_site"])),
            AvoidPattern("AAAAAA"),
            AvoidPattern("TTTTTT"),
            AvoidPattern("GGGGGG"),
            AvoidPattern("CCCCCC"),
            AvoidHairpins(),
            EnforceGCContent(mini=0.45, maxi=0.6, window=50)
        ],
        objectives=None
    )

    # SOLVE THE CONSTRAINTS
    problem.resolve_constraints()

    # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS
    print(problem.constraints_text_summary())

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)
    final_sequence = problem.sequence  # string

    return final_sequence


def valid_end_dna(five_prime_dna: str, three_prime_dna: str, enzyme: str,
                  fusion_sites_used_by_backbone: Tuple[str, ...]) -> None:
    class EndDna:
        def __init__(self, name: str, seq: str, the_frag_after_digest: int) -> None:
            self.name = name
            self.seq = Seq(seq)
            self.the_frag_after_digest = the_frag_after_digest

    end_dna_5: EndDna = EndDna("5'-sequence", five_prime_dna, -1)
    end_dna_3: EndDna = EndDna("3'-sequence", three_prime_dna, 0)

    biopy_enzyme = globals()[enzyme]
    print(f"Validating restriction enzyme sites in 5'-sequence and 3'-sequence:"
          f"\n 1. {enzyme} is found one time in 5'-sequence and 3'-sequence."
          f"\n 2. The fusion sites ,produced after digesting 5'-sequence and 3'-sequence using {enzyme} are also in the"
          f" specified fusion sites used by backbone.")
    validated_flag = True
    for end_dna in [end_dna_5, end_dna_3]:
        cuttable = True
        search_enzyme_sites: list = biopy_enzyme.search(end_dna.seq)
        if len(search_enzyme_sites) == 0:
            validated_flag = False
            cuttable = False
            warnings.warn(f"Warning! {enzyme} site is not found in {end_dna.name}.")
        elif len(search_enzyme_sites) > 1:
            validated_flag = False
            warnings.warn(f"Warning! {enzyme} site is found multiple times in {end_dna.name}.")

        # For example, BsaI.elucidate() outputs 'GGTCTCN^NNNN_N'
        fusion_site_length = len(re.sub('^.*\\^|_.*$', '', biopy_enzyme.elucidate()))
        if cuttable:
            if end_dna.name == "5'-sequence":
                fusion_site = biopy_enzyme.catalyze(end_dna.seq)[end_dna.the_frag_after_digest][:fusion_site_length]
            else:
                fusion_site = biopy_enzyme.catalyze(end_dna.seq)[end_dna.the_frag_after_digest][-fusion_site_length:]
            if fusion_site not in fusion_sites_used_by_backbone:
                validated_flag = False
                warnings.warn(f"Warning! The overhang produced by {enzyme} digestion of {end_dna.name}, "
                              f"i.e. {fusion_site}), is not found in the specified fusion sites used by backbone, "
                              f"i.e., {fusion_sites_used_by_backbone}.")
    if validated_flag:
        print("All validated! :) ")


def append_end_dna(fragment_with_fusion_sites: dict, five_prime_dna: str, three_prime_dna: str) -> dict:
    appended = copy.deepcopy(fragment_with_fusion_sites)
    for frag in appended.values():
        if frag["order"] == 0:
            # The parentheses specify the DNA additional to the DNA of the given input_draft amino acid sequence
            frag.update({"n_term_dna": "".join(["(", five_prime_dna, ")"])})
        if frag["order"] == len(appended) - 1:
            frag.update({"c_term_dna": "".join(["(", three_prime_dna, ")"])})

    return appended


def append_enzyme_sites_and_stuffer(fragment_with_fusion_sites: dict, enzyme: str, min_dna_frag_length: int) -> dict:
    biopy_enzyme = globals()[enzyme]
    dna_before_cutsite = re.sub("\\^.*$", "", biopy_enzyme.elucidate())
    appended = copy.deepcopy(fragment_with_fusion_sites)
    patterns = r"[()<>|]"

    for frag in appended.values():

        # For example, BsaI.elucidate() outputs 'GGTCTCN^NNNN_N', dna_to_add_nterm = 'GGTCTCN'
        dna_to_add_nterm = re.sub("N", "T", dna_before_cutsite)
        dna_to_add_nterm = "".join(["TTT", dna_to_add_nterm])

        # For example, BsaI.elucidate() outputs 'GGTCTCN^NNNN_N', dna_to_add_cterm = 'NGAGACC'
        dna_to_add_cterm = make_rev_compliment(dna_before_cutsite)
        dna_to_add_cterm = re.sub("N", "T", dna_to_add_cterm)
        dna_to_add_cterm = "".join([dna_to_add_cterm, "TTT"])

        # The parentheses specify the DNA additional to the DNA of the given input_draft amino acid sequence
        if frag["order"] == 0:
            # Do not need to append stuff to n_term_dna because this is taken care by append_end_dna()
            dna_to_add_nterm = ""
        if frag["order"] == len(appended) - 1:
            # Do not need to append stuff to c_term_dna because this is taken care by append_end_dna()
            dna_to_add_cterm = ""

        current_nterm_frag_dna_length = len(re.sub(patterns, "", frag.get("n_term_dna", ""))) + len(dna_to_add_nterm)
        current_cterm_frag_dna_length = len(re.sub(patterns, "", frag.get("c_term_dna", ""))) + len(dna_to_add_cterm)
        current_middle_frag_dna_length = len(frag.get("middle_aa")) * 3
        current_frag_dna_length = sum([current_nterm_frag_dna_length,
                                       current_cterm_frag_dna_length,
                                       current_middle_frag_dna_length])

        stuffer_dna_n = stuffer_dna_c = ""
        if current_frag_dna_length < min_dna_frag_length:
            stuffer_dna_n = STUFFER_REF[:(min_dna_frag_length - current_frag_dna_length) // 2]
            length_stuffer_c = (min_dna_frag_length - current_frag_dna_length - len(stuffer_dna_n))
            stuffer_dna_c = STUFFER_REF[len(stuffer_dna_n):len(stuffer_dna_n) + length_stuffer_c]

        frag.update({
            "n_term_dna": f"({stuffer_dna_n}{dna_to_add_nterm}){frag.get('n_term_dna', '')}",
            "c_term_dna": f"{frag.get('c_term_dna', '')}({dna_to_add_cterm}{stuffer_dna_c})",
        })

    return appended


def make_mutant_aa_fragments(fragment_n_and_c_term_dna: dict, mutations_0idx: list, linked_mutations_0idx: list,
                             codon_usage_table_path: str, positions_include_wt_aa_0idx: List[int]) -> list:
    patterns = r"[()<>|]"

    # Sort the partitioned sequence by order
    fragment_n_and_c_term_dna = dict(sorted(fragment_n_and_c_term_dna.items(), key=lambda x: x[1]['order']))
    fragments = [f for f in fragment_n_and_c_term_dna.keys()]
    start_positions = [sum([len(f) for f in fragments[0:i]]) for i in range(0, len(fragments) + 1)]

    mutations_0idx_reformat = {}
    if mutations_0idx:
        mutations_0idx_ = remove_linked_mutations_from_all_mutations(all_mutations=mutations_0idx,
                                                                     linked_mutations=linked_mutations_0idx)
        for d in mutations_0idx_:
            mut_new_format = []
            if d['position'] in positions_include_wt_aa_0idx:
                mut_new_format.append([(d['position'], 'wt')])
            for aa in d['aa']:
                mut_new_format.append([(d['position'], aa)])
            mutations_0idx_reformat.update({(d['position'],): mut_new_format})

    linked_mutations_reformat = {}
    if linked_mutations_0idx:
        for mut_set in linked_mutations_0idx:
            key = tuple([mut[1] for mut in mut_set])
            mut_set_format = [(mut[1], mut[2]) for mut in mut_set]
            include_wt = False
            if all([mut[1] in positions_include_wt_aa_0idx for mut in mut_set]):
                include_wt = True
            if key in linked_mutations_reformat:
                linked_mutations_reformat[key].append(mut_set_format)
            else:
                linked_mutations_reformat.update({key: [mut_set_format]})
                if include_wt:
                    linked_mutations_reformat[key].append([(mut[1], 'wt') for mut in mut_set])
        mutations_0idx_reformat.update(linked_mutations_reformat)

    # Initialize output_secret list
    output = []

    # Iterate over each fragment in the partitioned sequence
    for fragment, details in fragment_n_and_c_term_dna.items():
        fragment_start = start_positions[details['order']]
        fragment_end = start_positions[details['order'] + 1]

        # Initialize mut_positions_in_this_frag
        mut_positions_in_this_frag = []
        # Iterate over each mutation
        for positions, mutations in mutations_0idx_reformat.items():
            # Check if the mutation position is within the current fragment
            if fragment_start <= min(positions) and max(positions) < fragment_end:
                mut_positions_in_this_frag.append(mutations)

        for variant in product(*mut_positions_in_this_frag):
            # Generate the new fragment for mutation
            new_fragment_lst = list(fragment)

            # Generate the new middle aa for mutation
            search_start = re.search(details['middle_aa'], fragment)
            search_end = re.search(details['middle_aa'], fragment)
            middle_aa_start = search_start.start() if search_start else None
            middle_aa_end = search_end.end() if search_end else None
            new_middle_aas_lst = list(details['middle_aa'])
            new_n_term_dna = details['n_term_dna']
            new_c_term_dna = details['c_term_dna']

            mut_notations = []
            for mutation in variant:
                for position, aa in mutation:
                    if aa != 'wt':
                        # Add mutation notation
                        mut_notations.append(f"{fragment[position - fragment_start]}{position + 1}{aa}")

                        # Mutate the corresponding position on the fragment
                        new_fragment_lst[position - fragment_start] = aa

                        # Mutate the corresponding position on middle_aa
                        if middle_aa_start <= position - fragment_start < middle_aa_end:
                            new_middle_aas_lst[position - fragment_start - middle_aa_start] = aa

                        # position is not in middle aa but in the c-term dna
                        elif position - fragment_start >= middle_aa_end:
                            bases_before_fusion_site = re.sub('<.+$', '', details['c_term_dna'])
                            dna = re.sub(patterns, '', details['c_term_dna'])
                            codon_to_mutate_start = 3 * ((position - fragment_start) - middle_aa_end)
                            codon_to_mutate_end = 3 * ((position - fragment_start) - middle_aa_end + 1)
                            codon_to_mutate = dna[codon_to_mutate_start:codon_to_mutate_end]
                            if codon_to_mutate_end > len(bases_before_fusion_site):
                                # codon to mutate overlaps fusion site
                                overlap_bases_count = codon_to_mutate_end - len(bases_before_fusion_site)
                                overlap_codon_idxes = list(range(2, 2 - overlap_bases_count, -1))
                                fix_base = {
                                    idx: codon_to_mutate[idx] if idx in overlap_codon_idxes else None
                                    for idx in range(0, 3)
                                }
                                sel_codon = find_best_codon_by_usage(codon_usage_table_path, aa, fix_base)
                                if sel_codon != "Not found":
                                    new_c_term_dna = list(details['c_term_dna'])
                                    new_bases = [sel_codon[key] for key, v in fix_base.items() if v is None]
                                    new_c_term_dna[:len(new_bases)] = new_bases
                                    new_c_term_dna = ''.join(new_c_term_dna)
                                else:
                                    raise ValueError(
                                        f'Can not find codon for mutation {mutation} that fits the constraint of '
                                        f'fusion site!'
                                    )

                        # position is not in middle aa but in the n-term dna
                        elif position - fragment_start < middle_aa_start:
                            dna = re.sub(patterns, '', details['n_term_dna'])
                            bases_after_fusion_site = re.sub('^.+>', '', details['n_term_dna'])
                            codon_to_mutate_start = len(dna) - 3 * (middle_aa_start - (position - fragment_start))
                            codon_to_mutate_end = len(dna) - 3 * (middle_aa_start - (position - fragment_start) - 1)
                            codon_to_mutate = dna[codon_to_mutate_start:codon_to_mutate_end]
                            if len(dna) - codon_to_mutate_start > len(bases_after_fusion_site):
                                # codon to mutate overlaps fusion site
                                overlap_bases_count = len(dna) - codon_to_mutate_start - len(bases_after_fusion_site)
                                overlap_codon_idxes = list(range(0, overlap_bases_count, 1))
                                fix_base = {
                                    idx: codon_to_mutate[idx] if idx in overlap_codon_idxes else None
                                    for idx in range(0, 3)
                                }
                                sel_codon = find_best_codon_by_usage(codon_usage_table_path, aa, fix_base)
                                if sel_codon != "Not found":
                                    new_n_term_dna = list(details['n_term_dna'])
                                    new_bases = [sel_codon[key] for key, v in fix_base.items() if v is None]
                                    new_n_term_dna[-len(new_bases):] = new_bases
                                    new_n_term_dna = ''.join(new_n_term_dna)
                                else:
                                    raise ValueError(
                                        f'Can not find codon for mutation {mutation} that fits the constraint of '
                                        f'fusion site!'
                                    )

            new_fragment = ''.join(new_fragment_lst)
            new_middle_aas = ''.join(new_middle_aas_lst)

            if len(mut_notations) > 0:
                frag_name = '_'.join(sorted(list(set(mut_notations)), key=lambda x: re.sub('[A-Z]', '', x)))
            else:
                frag_name = "wild_type"
            frag_name = f"{fragment_start + 1}-{fragment_end}_{frag_name}"
            # Add the mutated fragment
            output.append({
                'name': frag_name,
                'fragment_start': fragment_start,
                'fragment_end': fragment_end,
                'fragment': new_fragment,
                'n_term_dna': new_n_term_dna,
                'c_term_dna': new_c_term_dna,
                'middle_aa': new_middle_aas
            })

    return output


def make_mutant_dna_fragments_from_mutant_aa_fragments(mutant_aa_fragments: list, codon_usage_table_path: str,
                                                       specify_dna_seq: str, enzyme: str) -> list:
    patterns = r"[()<>|]"
    # Note that the dict.copy() is shallow, if there is a nested list/etc in there changes will be applied to
    # both. IIRC. Deepcopy will avoid that.
    mutant_dna_fragments = copy.deepcopy(mutant_aa_fragments)
    output = []
    # specify_dna_seq needs to exactly encode the full aa sequence used for partitioning
    for fragment in mutant_dna_fragments:
        # make dna for middle_aa
        middle_aa_start = fragment['fragment_start'] + re.search(fragment['middle_aa'], fragment['fragment']).start()
        middle_aa_end = middle_aa_start + len(fragment['middle_aa'])
        fix_base = {0: None, 1: None, 2: None}
        mutant_middle_dna = "".join(
            [find_best_codon_by_usage(codon_usage_table_path, aa, fix_base) for aa in fragment['middle_aa']]
        )
        mutant_middle_dna = remove_enzyme_sites_in_dna(mutant_middle_dna, enzyme, codon_usage_table_path)
        mutant_middle_dna_ls = list(mutant_middle_dna)
        mutation_aa_indexes = []
        if len(specify_dna_seq) > 0:
            middle_dna = specify_dna_seq[middle_aa_start * 3: middle_aa_end * 3]
            wt_middle_aa = Seq(middle_dna).translate()
            mutation_aa_indexes = \
                [idx for idx, (aa_wt, aa_mut) in enumerate(zip(wt_middle_aa, fragment['middle_aa'])) if aa_wt != aa_mut]
            # replace wt codons with mutated codons
            middle_dna_ls = list(middle_dna)
            for mut_idx in mutation_aa_indexes:
                middle_dna_ls[mut_idx * 3:(mut_idx + 1) * 3] = mutant_middle_dna_ls[mut_idx * 3:(mut_idx + 1) * 3]
            middle_dna = "".join(middle_dna_ls)
        else:
            middle_dna = "".join(mutant_middle_dna_ls)
        fragment.update({'middle_dna': middle_dna})

        # concat n_term_dna, middle_dna, c_term_dna
        concatenated_dna = "".join([
            re.sub(patterns, "", fragment['n_term_dna']),
            re.sub(patterns, "", fragment['middle_dna']),
            re.sub(patterns, "", fragment['c_term_dna']),
        ])
        fragment.update({'concatenated_dna': concatenated_dna})

        validated = validate_mutant_dna_fragment_constraints(fragment, enzyme)
        if not validated:
            if len(specify_dna_seq) > 0:
                fragment = optimize_middle_dna_of_a_fragment_with_specify_dna_seq(
                    fragment, mutation_aa_indexes, enzyme, patterns, codon_usage_table_path
                )
        output.append(fragment)
    return output


def find_enzyme_sites_in_dna(dna_seq: str, enzyme: str, print_result: bool) -> List[tuple]:
    # Validate forward and reverse enzyme site only appear once each
    biopy_enzyme = globals()[enzyme]
    fw_enzyme_site = biopy_enzyme.site
    rv_enzyme_site = make_rev_compliment(fw_enzyme_site)
    locs = []
    for enzyme_site in [fw_enzyme_site, rv_enzyme_site]:
        enzyme_site_count = len(re.findall(enzyme_site, dna_seq))
        if enzyme_site_count != 0:
            locs_ = \
                [(m.start(0), m.end(0)) for m in re.finditer(enzyme_site, dna_seq)]
            locs = locs + locs_
            if print_result:
                print(f' {enzyme} site {enzyme_site} is found {enzyme_site_count} times!')
                print(f'Locations: {locs_}')
    return locs


def remove_enzyme_sites_in_dna(dna_seq: str, enzyme: str, codon_usage_table_path: str) -> str:
    if len(dna_seq) % 3 != 0:
        raise ValueError('Partial codon. DNA sequence is not a multiple of three.')

    # assume first nt as the translation start
    locs = find_enzyme_sites_in_dna(dna_seq, enzyme, False)
    codon_usage: dict = generate_aa2codon_dict(codon_usage_table_path)
    sel_dna_seq = dna_seq
    sel_dna_seq_enzyme_locs = locs
    for loc in locs:
        optimised = False
        aa_idxes_of_codons_to_be_changed = list(range(loc[0] // 3, loc[1] - 1 // 3 + 1))
        for aa_idx in aa_idxes_of_codons_to_be_changed:
            if optimised: break
            aa = Seq(dna_seq)[aa_idx * 3: (aa_idx + 1) * 3].translate()
            codons_sorted_by_usage = [cod for cod, frac in
                                      sorted(zip(codon_usage[aa]['codon'], codon_usage[aa]['relative_frequency']),
                                             key=lambda pair: pair[0])]

            # Try each codon until complying constraints
            for new_codon in codons_sorted_by_usage:
                new_dna_seq = MutableSeq(sel_dna_seq)
                new_dna_seq[aa_idx * 3: (aa_idx + 1) * 3] = new_codon
                new_enzyme_locs = find_enzyme_sites_in_dna(str(new_dna_seq), enzyme, False)
                if len(new_enzyme_locs) < len(sel_dna_seq_enzyme_locs):
                    optimised = True
                    sel_dna_seq = new_dna_seq
                    sel_dna_seq_enzyme_locs = new_enzyme_locs
                    break

        if not optimised:
            print(f'Warning! {enzyme} site at {loc} can not be removed!')
            sel_dna_seq_enzyme_locs = find_enzyme_sites_in_dna(str(sel_dna_seq), enzyme, False)
            print(f'The optimized sequence has {len(sel_dna_seq_enzyme_locs)} {enzyme} sites '
                  f'at these locations: {sel_dna_seq_enzyme_locs}')

    return str(sel_dna_seq)


def optimize_middle_dna_of_a_fragment_with_specify_dna_seq(dna_fragment: dict, mutation_aa_indexes: list, enzyme: str,
                                                           patterns: str, codon_usage_table_path: str) -> dict:
    dna_fragment_working_copy = copy.deepcopy(dna_fragment)
    # Validate dna sequence constraints
    codon_usage: dict = generate_aa2codon_dict(codon_usage_table_path)

    # Try each codon for each mutation position until complying constraints
    passed = False
    if len(mutation_aa_indexes) == 0:
        if validate_mutant_dna_fragment_constraints(dna_fragment_working_copy, enzyme):
            return dna_fragment_working_copy
        else:
            warnings.warn(f'found {enzyme} site in concatenated wt dna fragment {dna_fragment_working_copy["name"]}')
    else:
        for mut_idx in mutation_aa_indexes:
            aa = dna_fragment_working_copy['middle_aa'][mut_idx]
            codons_sorted_by_usage = [cod for cod, frac in
                                      sorted(zip(codon_usage[aa]['codon'], codon_usage[aa]['relative_frequency']),
                                             key=lambda pair: pair[0])]

            # Try each codon until complying constraints
            for new_codon in codons_sorted_by_usage:
                middle_dna_ls = list(dna_fragment_working_copy['middle_dna'])
                middle_dna_ls[mut_idx * 3:(mut_idx + 1) * 3] = new_codon
                dna_fragment_working_copy.update({'middle_dna': "".join(middle_dna_ls)})

                # concat n_term_dna, middle_dna, c_term_dna
                concatenated_dna = "".join([
                    re.sub(patterns, "", dna_fragment_working_copy['n_term_dna']),
                    re.sub(patterns, "", dna_fragment_working_copy['middle_dna']),
                    re.sub(patterns, "", dna_fragment_working_copy['c_term_dna']),
                ])
                dna_fragment_working_copy.update({'concatenated_dna': concatenated_dna})
                if validate_mutant_dna_fragment_constraints(dna_fragment_working_copy, enzyme):
                    passed = True
                    break
            if passed:
                return dna_fragment_working_copy
    if not passed:
        msg = f"{dna_fragment_working_copy['name']} Unable to find a dna sequence that complies the constraints!"
        warnings.warn(msg)
        return dna_fragment


def validate_mutant_dna_fragment_constraints(mutant_dna_fragment: dict, enzyme: str) -> bool:
    # Sequences that cause problems for gBlocks Gene Fragment synthesis are those having extremely low or high GC
    # content (less than 25% and greater than 75%), homopolymeric runs of 10 or more As and Ts or 6 or more Gs and
    # Cs. Other structural motifs such as repeats or hairpins may also influence our ability to construct gBlocks
    # Gene Fragments. The decision to accept or reject an order for a particular gBlocks Gene Fragment sequence is
    # complex and multifactorial. (source: IDT FAQs)

    passed = validate_dna_seq_enzyme_constraints(dna_seq=mutant_dna_fragment['concatenated_dna'], enzyme=enzyme)
    if not passed:
        print(f'Validation failed for fragment {mutant_dna_fragment["name"]}...')
        print(f'Forward and reverse {enzyme} site should only appear once each.')

    return passed


def validate_dna_seq_enzyme_constraints(dna_seq: str, enzyme: str) -> bool:
    # Validate forward and reverse enzyme site only appear once each
    passed = True
    biopy_enzyme = globals()[enzyme]
    fw_enzyme_site = biopy_enzyme.site
    rv_enzyme_site = make_rev_compliment(fw_enzyme_site)
    fw_enzyme_site_count = len(re.findall(fw_enzyme_site, dna_seq))
    rv_enzyme_site_count = len(re.findall(rv_enzyme_site, dna_seq))
    if fw_enzyme_site_count != 1:
        passed = False
        locations = \
            [(m.start(0), m.end(0)) for m in re.finditer(fw_enzyme_site, dna_seq)]
        print(f'Forward {enzyme} site {fw_enzyme_site} is found {fw_enzyme_site_count} times!')
        print(f'Locations: {locations}')
    if rv_enzyme_site_count != 1:
        passed = False
        locations = \
            [(m.start(0), m.end(0)) for m in re.finditer(rv_enzyme_site, dna_seq)]
        print(f'Reverse {enzyme} site {rv_enzyme_site} is found {rv_enzyme_site_count} times!')
        print(f'Locations: {locations}')

    return passed


def validate_mutant_dna_fragments_constraints(mutant_dna_fragments: list, enzyme: str) -> bool:
    all_passed = True
    for mutant_dna_fragment in mutant_dna_fragments:
        if not validate_mutant_dna_fragment_constraints(mutant_dna_fragment, enzyme):
            all_passed = False

    return all_passed


def export_module_ordering_sheet(gene_abbreviation: str, mutant_dna_fragments: list, row_range: Tuple[str, str],
                                 column_range: Tuple[int, int], enzyme: str, output_file: str) -> None:
    print('Exporting DNA module ordering sheet')
    well_positions = multi_well_plate_position_generator(row_range, column_range)
    mutant_dna_fragments_working_copy = copy.deepcopy(mutant_dna_fragments)
    last_frag_identity = ''
    frag_id_count = 1
    for idx, frag in enumerate(mutant_dna_fragments_working_copy):
        enzyme_recognition_pat = re.sub('N','',ENZYME_INFO[enzyme]['recognition_site'])
        enzyme_recognition_pat_rc = str(Seq(enzyme_recognition_pat).reverse_complement())
        labeled_dna = re.sub(enzyme_recognition_pat, f"<{enzyme_recognition_pat}>", frag['concatenated_dna'])
        labeled_dna = re.sub(enzyme_recognition_pat_rc, f"<{enzyme_recognition_pat_rc}>", labeled_dna)
        frag_identity = re.sub("_wild_type|_[A-Z][0-9]+[A-Z]", "", frag["name"])
        if frag_identity == last_frag_identity:
            frag_id_count += 1
        else:
            frag_id_count = 1
        last_frag_identity = frag_identity
        short_name = f'{gene_abbreviation}_{re.sub("_wild_type|_[A-Z][0-9]+[A-Z]", "", frag["name"])}_{frag_id_count}'
        frag.update({'Well Position': well_positions[idx % len(well_positions)],
                     'Labeled Sequence': labeled_dna,
                     'Short Name': short_name})

    with pd.ExcelWriter(output_file) as writer:
        for count in range(0, len(mutant_dna_fragments_working_copy) // len(well_positions) + 1):
            plate = mutant_dna_fragments_working_copy[len(well_positions) * count:len(well_positions) * (count + 1)]
            df = pd.DataFrame(plate)[['Well Position',
                                      'name',
                                      'concatenated_dna',
                                      'Labeled Sequence',
                                      'fragment',
                                      'Short Name']]
            df.rename(columns={'fragment': 'Amino Acid Sequence',
                               'name': 'Sequence Name',
                               'concatenated_dna': 'Sequence'},
                      inplace=True)
            df.to_excel(writer, sheet_name=f'plate{count + 1}', index=False)
        print('--------------------------------------------------------------------------------------------------')
        print(f'\033[1m DNA module ordering sheet is exported at:\n {output_file}\033[0m')


def import_mutant_dna_fragments_from_module_sheet(module_sheet_path: str) -> List[dict]:
    plate_map_df = load_module_sheet(module_sheet_path)
    plate_map_df.rename(columns={'Amino Acid Sequence': 'fragment',
                                 'Sequence Name': 'name',
                                 'Sequence': 'concatenated_dna'},
                        inplace=True)
    imported_mutant_dna_fragments = plate_map_df.to_dict(orient='records')
    return imported_mutant_dna_fragments


def process_fragment_with_fusion_sites(input_params: dict, fragment_with_fusion_sites: dict,
                                       check_seq_complexity_idt: bool, product_type: str,
                                       supress_output: bool, min_dna_frag_length: int,
                                       positions_include_wt_aa_0idx: List[int]) -> tuple:
    all_mutations_0idx, linked_mutations_0idx = prepare_0idx_mutations(input_params['mutations_1idx'],
                                                                       input_params['linked_mutations_1idx'])

    if not supress_output:
        print('Processing the partitioned fragments...')
        print("  Concatenating user-specified 5' and 3' dna sequence...")
    fragment_with_three_and_five_prime_dna = append_end_dna(
        fragment_with_fusion_sites=fragment_with_fusion_sites,
        five_prime_dna=input_params['five_prime_dna'],
        three_prime_dna=input_params['three_prime_dna'],
    )
    if not supress_output:
        print(f"  Concatenating {input_params['enzyme']} sites and stuffers...")
    fragment_n_and_c_term_dna = append_enzyme_sites_and_stuffer(
        fragment_with_fusion_sites=fragment_with_three_and_five_prime_dna,
        enzyme=input_params['enzyme'],
        min_dna_frag_length=min_dna_frag_length
    )
    if not supress_output:
        print(f"  Generating individual mutant aa fragments...")
    mutant_aa_fragments = make_mutant_aa_fragments(
        fragment_n_and_c_term_dna=fragment_n_and_c_term_dna,
        mutations_0idx=all_mutations_0idx,
        linked_mutations_0idx=linked_mutations_0idx,
        codon_usage_table_path=input_params['codon_usage_table_path'],
        positions_include_wt_aa_0idx=positions_include_wt_aa_0idx
    )
    if not supress_output:
        print(f"  Translating mutant aa fragments to dna fragments...")
    mutant_dna_fragments = make_mutant_dna_fragments_from_mutant_aa_fragments(
        mutant_aa_fragments=mutant_aa_fragments,
        codon_usage_table_path=input_params['codon_usage_table_path'],
        specify_dna_seq=input_params['fix_wt_dna_sequence'],
        enzyme=input_params['enzyme']
    )
    if not supress_output:
        print("  Validating dna fragment constraints (enzyme sites)...")
    if not validate_mutant_dna_fragments_constraints(mutant_dna_fragments,
                                                     input_params['enzyme']) and not supress_output:
        warnings.warn("   Mutant dna fragments failed constraint validation!!!!!!")
    else:
        if not supress_output:
            print("   Passed validation of mutant dna fragments constraints!")

    dna = [{f['name']: f['concatenated_dna']} for f in mutant_dna_fragments]
    if check_seq_complexity_idt:
        passed = validate_seq_complexity_idt(
            dna=dna, idt_credentials=IdtCredentials, product_type=product_type
        )
        if not passed and not supress_output:
            warnings.warn(f"   Mutant dna fragments failed IDT complexity check for {product_type} !!!!!!")
        else:
            if not supress_output:
                print(f"   Passed IDT complexity check for {product_type}!!!!!!")

    return mutant_aa_fragments, mutant_dna_fragments


def validate_partitioned_fragments_by_insilico_assembly(mutant_dna_fragments: List[dict], sample_number: int,
                                                        wt_seq: str, enzyme: str, five_prime_dna: str,
                                                        three_prime_dna: str, coding_start: Union[int, None]) -> bool:
    print('Validating mutant dna fragments by in-silico assembly...')
    validated = True
    count = 0
    while count < sample_number:
        random_combi = random_fragment_combinations(mutant_dna_fragments=mutant_dna_fragments)
        assembled_frag = assemble_fragments(a_set_of_dna_fragments=random_combi, enzyme=enzyme, wt_seq=wt_seq,
                                            five_prime_dna=five_prime_dna, three_prime_dna=three_prime_dna)
        validated = validated & assembled_frag.validate_assembled_fragment()
        print(f'Random combi in-silico assembly: {assembled_frag.name}, \n  Validation Passed: {validated}')
        count += 1
    return validated


def post_partition_processing(input_file_path: str, best_partitions_by_cut_number_file: str, cut_number: int,
                              positions_include_wt_aa_0idx: List[int], check_seq_complexity_idt: bool,
                              product_type: str, validate_sample_number: int, validate_coding_start: Union[None, int],
                              min_dna_frag_length: int, cost_per_nt: float) -> tuple:
    print('\033[1m===========================================================================================\033[0m')
    print(f'                               \033[1m GENERATING MODULES \033[0m '
          f'\n                               Number of fragments: {cut_number + 1} ')
    print('\033[1m===========================================================================================\033[0m')

    input_params = load_input_params(input_file_path=input_file_path, supress_output=True)

    with open(best_partitions_by_cut_number_file) as f:
        best_partitions_by_cut_number = json.load(f)

    sel_res = None
    for res_per_cut in best_partitions_by_cut_number["best_partitions_by_cut_number"]:
        if res_per_cut["number_of_cuts"] == cut_number:
            sel_res = res_per_cut
    if not sel_res:
        raise ValueError(f'Unable to find the results of the desired cut number: {cut_number} in the specified '
                         f'best_partitions_by_cut_number_file: {best_partitions_by_cut_number_file}')
    if len(sel_res["sel_partitions"]) == 0:
        raise ValueError(f'No partitions are found for the selected cut number:{cut_number}!')

    for partition in sel_res["sel_partitions"]:
        print(f'\n Processing partition:\n'
              f' partition: {partition["partition"]}\n'
              f' ligation_fidelity: {partition["ligation_fidelity"]}\n'
              f' fragment_length_unevenness: {partition["fragment_length_unevenness"]}\n'
              f' fusion_sites: {partition["fusion_sites"]}\n'
              f' cost: {partition["cost"]}\n'
              f' fragments: \n {partition["fragments"]}\n')

        fragment_with_fusion_sites = partition["fragment_with_fusion_sites"]

        mutant_aa_fragments, mutant_dna_fragments = \
            process_fragment_with_fusion_sites(input_params, fragment_with_fusion_sites, check_seq_complexity_idt,
                                               product_type, supress_output=False,
                                               min_dna_frag_length=min_dna_frag_length,
                                               positions_include_wt_aa_0idx=positions_include_wt_aa_0idx)

        if validate_partitioned_fragments_by_insilico_assembly(mutant_dna_fragments=mutant_dna_fragments,
                                                               sample_number=validate_sample_number,
                                                               wt_seq=input_params['s'],
                                                               enzyme=input_params['enzyme'],
                                                               five_prime_dna=input_params['five_prime_dna'],
                                                               three_prime_dna=input_params['three_prime_dna'],
                                                               coding_start=validate_coding_start):
            mtp_format = input_params['module_plate_format']
            if mtp_format not in PLATE_FORMATS.keys():
                raise ValueError(f'Invalid MTP format. MTP format muts be one of these: {list(PLATE_FORMATS.keys())}')

            outfile_path = path.join(path.dirname(path.dirname(best_partitions_by_cut_number_file)), 'results',
                                     f'order_modules_{cut_number+1}fragments.xlsx')
            export_module_ordering_sheet(
                gene_abbreviation=input_params['gene_name'],
                mutant_dna_fragments=mutant_dna_fragments,
                row_range=PLATE_FORMATS[mtp_format]['rows'],
                column_range=PLATE_FORMATS[mtp_format]['columns'],
                enzyme=input_params['enzyme'],
                output_file=outfile_path
            )

            exact_cost = find_exact_cost(mutant_dna_fragments=mutant_dna_fragments,
                                         price_per_base=cost_per_nt)
            print(
                f'\n\033[1m'
                f'\n================================================================================================'
                f'\n                                SUCCESSFULLY GENERATED MODULES!                       '
                f'\n                                Number of fragments: {cut_number+1}                       '
                f'\n                                 Estimated cost: {exact_cost} €'
                f'\n================================================================================================'
                f'\033[0m'
            )

            return mutant_aa_fragments, mutant_dna_fragments, outfile_path

        else:
            raise ValueError('validate_partitioned_fragments_by_insilico_assembly() failed')
