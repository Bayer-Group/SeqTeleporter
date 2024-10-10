import ast
from os import environ, listdir
from os.path import join, dirname, abspath, basename
from dotenv import load_dotenv
from Bio.Seq import Seq
import pandas as pd
import re
import math
from typing import Tuple
from python_codon_tables.python_codon_tables import _tables_dir as codon_usage_tbl_dir

from seqteleporter.utils.utils import (include_linked_mutations_into_mutations, annotate_sequence_mutations,
                                       prepare_0idx_mutations)


def transform_input_excel_sheet_to_text_input(input_table_path: str) -> Tuple[str, dict]:

    input_param_specs = \
        pd.read_excel(input_table_path, sheet_name='parameter_guide', header=0, index_col=None)\
        .to_dict(orient='records')
    input_param_specs_dict = {i['parameter']: i for i in input_param_specs}
    with pd.option_context('future.no_silent_downcasting', True):
        inputs_dict = (
            pd.read_excel(input_table_path, sheet_name='input_parameters', header=None, index_col=0, usecols=[0, 1])
            .transpose()
            .fillna({'FIX_DNA_SEQUENCE': '',
                     'FUSION_SITES_USED_BY_BACKBONE': '()',
                     'DNA_5_PRIME': '',
                     'DNA_3_PRIME': '',
                     'ALLOWED_CUT_POSITIONS': '[]',
                     'BACKBONE_SIZE': 3000
                     })
            .transpose()
            .to_dict()[1]
        )
    for k, v in inputs_dict.items():
        if type(v) == str: inputs_dict.update({k: re.sub(' ', '', v)})

    for k, v in inputs_dict.items():
        if input_param_specs_dict[k]['required'] == 'yes' and type(v) in [float, int] and math.isnan(v):
            raise ValueError(f'Error in input file! Input parameter {k} should not be empty!')
        if input_param_specs_dict[k]['required'] == 'yes' and v == '':
            raise ValueError(f'Error in input file! Input parameter {k} should not be empty!')

    if inputs_dict['FUSION_SITES_USED_BY_BACKBONE'] != '()':
        inputs_dict.update(
            {'FUSION_SITES_USED_BY_BACKBONE': tuple(inputs_dict['FUSION_SITES_USED_BY_BACKBONE'].split(','))})

    input_mutation_table = pd.read_excel(input_table_path, sheet_name='input_mutations', header=0, index_col=None)
    for f in listdir(join(dirname(dirname(__file__)), 'data', 'neb_fidelity_data')):
        if re.search(inputs_dict['FIDELITY_DATA'], f):
            inputs_dict.update({'FIDELITY_DATA': f})

    check_that_the_wt_aa_in_input_muts_are_correct(input_mutation_table=input_mutation_table,
                                                   wt_seq=inputs_dict['SEQUENCE'])

    linked_mut_sets_1dx = []
    for idx, df in input_mutation_table[~input_mutation_table.linked.isna()].groupby(['linked']):
        df['position_0idx'] = df['position'] - 1
        df['mut_tup'] = df.apply(lambda row: (inputs_dict['SEQUENCE'][row['position_0idx']],
                                              row['position'],
                                              row['mut_aa']),
                                 axis=1)
        linked_muts_1dx = tuple([mut_tup for mut_tup in df['mut_tup']])
        linked_mut_sets_1dx.append(linked_muts_1dx)

    mutations = input_mutation_table[['position', 'mut_aa']].rename(columns={'mut_aa': 'aa'})\
        .groupby(['position'])\
        .agg({'aa': lambda x: list(x.drop_duplicates())})\
        .reset_index()\
        .to_dict('records')

    inputs_dict.update({'MUTATIONS': mutations, 'LINKED_MUTATIONS': linked_mut_sets_1dx})

    out_file_path = join(dirname(input_table_path), f'{re.sub("xlsx","txt",basename(input_table_path))}')
    with open(out_file_path, 'w') as f:
        for k, v in inputs_dict.items():
            f.write(f'{k}={v}\n')

    return out_file_path, inputs_dict


def check_that_the_wt_aa_in_input_muts_are_correct(input_mutation_table: pd.DataFrame, wt_seq: str):
    for idx_, row in input_mutation_table[['position', 'wt_aa']].drop_duplicates().iterrows():
        if wt_seq[int(row['position'])-1] != row['wt_aa']:
            raise ValueError(f"Input mutation validation failed! "
                             f"Position {row['position']} has inconsistent wt amino acid in "
                             f"input_parameters: {wt_seq[int(row['position'])-1]} and input_mutations: {row['wt_aa']}")


def load_input_params(input_file_path: str, supress_output: bool) -> dict:
    load_dotenv(input_file_path, override=True)
    gene_name = environ.get("GENE_NAME")
    allowed_cut_positions_1idx = ast.literal_eval(environ.get("ALLOWED_CUT_POSITIONS"))
    five_prime_dna = environ.get("DNA_5_PRIME")
    three_prime_dna = environ.get("DNA_3_PRIME")
    fix_wt_dna_sequence = environ.get("FIX_DNA_SEQUENCE")
    fusion_sites_used_by_backbone = ast.literal_eval(environ.get("FUSION_SITES_USED_BY_BACKBONE"))
    host = environ.get("HOST")
    fidelity_data = environ.get("FIDELITY_DATA")
    fidelity_data_path = join(dirname(dirname(abspath(__file__))), 'data', 'neb_fidelity_data', fidelity_data)
    sequence = environ.get("SEQUENCE")
    mutations_1idx = ast.literal_eval(environ.get("MUTATIONS"))
    linked_mutations_1idx = ast.literal_eval(environ.get("LINKED_MUTATIONS"))
    cut_number_range = ast.literal_eval(environ.get("CUT_NUMBER_RANGE"))
    min_aa_length = int(environ.get("MIN_AA_LENGTH"))
    max_cost = float(environ.get("MAX_COST"))
    max_length_unevenness = float(environ.get("MAX_LENGTH_UNEVENNESS"))
    min_ligation_fidelity = float(environ.get("MIN_LIGATION_FIDELITY"))
    satisfaction_ligation_fidelity = float(environ.get("SATISFACTION_LIGATION_FIDELITY"))
    cost_per_nt = float(environ.get("COST_PER_BP"))
    provider_min_frag_len = int(environ.get("MIN_DNA_LENGTH"))
    provider_max_frag_len = int(environ.get("MAX_DNA_LENGTH"))
    module_plate_format = int(environ.get("MODULE_PLATE_FORMAT"))
    enzyme = environ.get("ENZYME")

    if not supress_output:
        print(f'GENE_NAME={gene_name}')
        print(f'ALLOWED_CUT_POSITIONS={allowed_cut_positions_1idx}')
        print(f'DNA_5_PRIME={five_prime_dna}')
        print(f'DNA_3_PRIME={three_prime_dna}')
        print(f'FIX_DNA_SEQUENCE={fix_wt_dna_sequence}')
        print(f'FUSION_SITES_USED_BY_BACKBONE={fusion_sites_used_by_backbone}')
        print(f'HOST={host}')
        print(f'FIDELITY_DATA_PATH={fidelity_data_path}')
        print(f'SEQUENCE={sequence}')
        print(f'MUTATIONS={mutations_1idx}')
        print(f'LINKED_MUTATIONS={linked_mutations_1idx}')
        print(f'CUT_NUMBER_RANGE={cut_number_range}')
        print(f'MIN_AA_LENGTH={min_aa_length}')
        print(f'MAX_COST={max_cost}')
        print(f'COST_PER_BP={cost_per_nt}')
        print(f'MIN_DNA_LENGTH={provider_min_frag_len}')
        print(f'MAX_DNA_LENGTH={provider_max_frag_len}')
        print(f'MAX_LENGTH_UNEVENNESS={max_length_unevenness}')
        print(f'MIN_LIGATION_FIDELITY={min_ligation_fidelity}')
        print(f'SATISFACTION_LIGATION_FIDELITY={satisfaction_ligation_fidelity}')
        print(f'MODULE_PLATE_FORMAT={module_plate_format}')
        print(f'ENZYME={enzyme}')

    # validate inputs
    codon_usage_table_path = ""
    for f in listdir(codon_usage_tbl_dir):
        if re.search(host, f):
            codon_usage_table_path = join(codon_usage_tbl_dir, f)
            break
    if codon_usage_table_path == "":
        raise ValueError(f'Unable to find a codon usage table for the provided host {host}.\n'
                         f'Here are the available codon usage data in the codon usage data folder '
                         f'{codon_usage_tbl_dir}:\n'
                         f'{[f for f in listdir(codon_usage_tbl_dir) if re.match(".*csv$", f)]}')

    params = dict(
        s=sequence,
        mutations_1idx=mutations_1idx,
        linked_mutations_1idx=linked_mutations_1idx,
        cut_number_range=cut_number_range,
        fidelity_data_path=fidelity_data_path,
        min_aa_length=min_aa_length,
        fusion_sites_used_by_backbone=fusion_sites_used_by_backbone,
        max_cost=max_cost,
        max_unevenness=max_length_unevenness,
        min_ligation_fidelity=min_ligation_fidelity,
        satisfaction_fidelity=satisfaction_ligation_fidelity,
        codon_usage_table_path=codon_usage_table_path,
        host=host,
        allowed_cut_positions_1idx=allowed_cut_positions_1idx,
        five_prime_dna=five_prime_dna,
        three_prime_dna=three_prime_dna,
        fix_wt_dna_sequence=fix_wt_dna_sequence,
        cost_per_nt=cost_per_nt,
        provider_min_frag_len=provider_min_frag_len,
        provider_max_frag_len=provider_max_frag_len,
        gene_name=gene_name,
        module_plate_format=module_plate_format,
        enzyme=enzyme
    )

    return params


def validate_input_params(input_file_path: str):
    print('\033[1mValidating your input file...!\033[0m')
    input_params = load_input_params(input_file_path=input_file_path, supress_output=True)
    validate_specified_dna(input_params)
    check_that_no_wt_aa_is_misassigned_as_mutation(input_params)
    print('\033[1mInput file passed validation!\033[0m')


def validate_specified_dna(input_params: dict):

    # check that the specified dna seq indeed encodes the given aa seq
    if len(input_params['fix_wt_dna_sequence']) > 0:
        if str(Seq(input_params['fix_wt_dna_sequence']).translate()) != input_params['s']:
            raise ValueError('The specified dna seq (FIX_DNA_SEQUENCE) does not encode the aa seq (SEQUENCE)!')


def check_that_no_wt_aa_is_misassigned_as_mutation(input_params: dict):
    all_mutations_1idx = include_linked_mutations_into_mutations(input_params['mutations_1idx'],
                                                                 input_params["linked_mutations_1idx"])
    caught_mis_assigned_wt_aas = []
    for mut in all_mutations_1idx:
        unique_mut_aas = list(set(mut['aa']))
        wt_aa = input_params['s'][mut['position']-1]
        if wt_aa in unique_mut_aas:
            caught_mis_assigned_wt_aas.append({'position': mut['position'], 'wt_aa': wt_aa})
    if len(caught_mis_assigned_wt_aas) != 0:
        raise ValueError(f'Found incorrect mutations: {caught_mis_assigned_wt_aas} in the input file.'
                         f'Please correct the mutations and try again.')


def show_input_seq_info(input_file_path):
    input_params = load_input_params(input_file_path=input_file_path, supress_output=True)
    all_mutations_0idx, linked_mutations_0idx = prepare_0idx_mutations(
        input_params['mutations_1idx'], input_params['linked_mutations_1idx']
    )
    s_annotated = annotate_sequence_mutations(input_params['s'], all_mutations_0idx, linked_mutations_0idx)
    print(f'\033[1m ------------------- \033[0m')
    print(f'\033[1m Here is your input: \033[0m')
    print(f'\033[1m ------------------- \033[0m')
    print(f'\033[1mGENE_NAME = \033[0m {input_params["gene_name"]}')
    print(f'\033[1mSequence (red: independent mutation positions; other colors: linked mutation positions) = \033[0m'
          f'\n{"".join(s_annotated.values())}')
    print(f'\033[1mSequence Length = \033[0m {len(input_params["s"])}')
    print(f'\033[1mMutations (1-indexed) = \033[0m {input_params["mutations_1idx"]}')
    print(f'\033[1mLinked mutations (1-indexed) = \033[0m {input_params["linked_mutations_1idx"]}')
    print(f'\033[1mALLOWED_CUT_POSITIONS=\033[0m{input_params["allowed_cut_positions_1idx"]}')
    print(f'\033[1mDNA_5_PRIME=\033[0m{input_params["five_prime_dna"]}')
    print(f'\033[1mDNA_3_PRIME=\033[0m{input_params["three_prime_dna"]}')
    print(f'\033[1mFIX_DNA_SEQUENCE=\033[0m{input_params["fix_wt_dna_sequence"]}')
    print(f'\033[1mFUSION_SITES_USED_BY_BACKBONE=\033[0m{input_params["fusion_sites_used_by_backbone"]}')
    print(f'\033[1mHOST=\033[0m{input_params["host"]}')
    print(f'\033[1mFIDELITY_DATA_PATH=\033[0m{input_params["fidelity_data_path"]}')
    print(f'\033[1mMUTATIONS=\033[0m{input_params["mutations_1idx"]}')
    print(f'\033[1mLINKED_MUTATIONS=\033[0m{input_params["linked_mutations_1idx"]}')
    print(f'\033[1mCUT_NUMBER_RANGE=\033[0m{input_params["cut_number_range"]}')
    print(f'\033[1mMIN_AA_LENGTH=\033[0m{input_params["min_aa_length"]}')
    print(f'\033[1mMAX_COST=\033[0m{input_params["max_cost"]}')
    print(f'\033[1mCOST_PER_BP=\033[0m{input_params["cost_per_nt"]}')
    print(f'\033[1mMIN_DNA_LENGTH=\033[0m{input_params["provider_min_frag_len"]}')
    print(f'\033[1mMAX_DNA_LENGTH=\033[0m{input_params["provider_max_frag_len"]}')
    print(f'\033[1mMAX_LENGTH_UNEVENNESS=\033[0m{input_params["max_unevenness"]}')
    print(f'\033[1mMIN_LIGATION_FIDELITY=\033[0m{input_params["min_ligation_fidelity"]}')
    print(f'\033[1mSATISFACTION_LIGATION_FIDELITY=\033[0m{input_params["satisfaction_fidelity"]}')
    print(f'\033[1mMODULE_PLATE_FORMAT=\033[0m{input_params["module_plate_format"]}')
    print(f'\033[1mENZYME=\033[0m{input_params["enzyme"]}')




