from datetime import date
import json
from os import path
import pandas as pd
import shutil

from seqteleporter.utils.load_input_params import (validate_input_params, transform_input_excel_sheet_to_text_input,
                                                   show_input_seq_info)
from seqteleporter.post_partition_processor.post_partition_processor import post_partition_processing
from seqteleporter.partitioner.compute_best_partitions import (compute_best_partitions,
                                                               prepare_compute_best_partitions_params)
from seqteleporter.fragment_assembler.plate_mapper import make_and_validate_plate_mapping_sheet, load_module_sheet
from seqteleporter.fragment_assembler.fragment_assembler import generate_all_possible_variants_from_modules


"""
This script contains functions to generate and optimize modules for a sequencing project, as well as to assemble those modules and generate robot instructions for further processing.

The main functionalities include:
1. `generate_and_optimize_ready_to_click_modules`: 
   - Takes an input Excel file containing sequencing data, processes it to validate parameters, and computes the best partitions based on specified cut numbers.
   - Generates mutant amino acid and DNA fragments, saves them in JSON format, and copies the original input file to the output directory for reference.
   - Returns a list of output file paths for the generated modules.

2. `assemble_modules_and_generate_robot_instruction`: 
   - Reads desired variants from an input Excel file and processes them to create a plate mapping sheet.
   - Generates possible variants based on module names or specified mutations and names.
   - Validates the generated plate mapping sheet against the provided parameters such as plate format and plasmid IDs.

Dependencies:
- `datetime`: For handling date-related functionalities.
- `json`: For saving data in JSON format.
- `os.path`: For file path manipulations.
- `pandas`: For data handling and processing.
- `shutil`: For file operations such as copying files.
- Custom utility functions from `seqteleporter` package for input validation, partition computation, and fragment assembly.

Usage:
- Call `generate_and_optimize_ready_to_click_modules` with the path to the input Excel file to generate modules.
- Call `assemble_modules_and_generate_robot_instruction` with paths to the input Excel file and the generated modules to create a mapping sheet for robotic processing.

Note: Ensure all necessary input files are in the correct format as expected by the utility functions to avoid runtime errors.
"""


def generate_and_optimize_ready_to_click_modules(input_table_path: str) -> list[str]:
    """
    Generates and optimizes modules for a sequencing project based on the input Excel file.

    This function processes the input data to validate parameters, computes the best partitions
    for the specified cut numbers, and generates mutant amino acid and DNA fragments. The results
    are saved in JSON format, and the original input file is copied to the output directory for reference.

    Args:
        input_table_path (str): The file path to the input Excel file containing sequencing data.

    Returns:
        list: A list of output file paths for the generated modules.
    """
    print('\033[1m=============================================================================================\033[0m')
    print('                               \033[1m RUN STARTED! \033[0m ')
    print('\033[1m=============================================================================================\033[0m')
    text_input_path, inputs_dict_ = transform_input_excel_sheet_to_text_input(input_table_path)
    validate_input_params(text_input_path)
    show_input_seq_info(text_input_path)
    params = prepare_compute_best_partitions_params(text_input_path)
    best_partitions_by_cut_number, best_partitions_by_cut_number_file_path = compute_best_partitions(**params)
    outfile_paths = []
    for cut_number in range(params['cut_number_range'][0], params['cut_number_range'][1]):
        mutant_aa_fragments, mutant_dna_fragments, outfile_path = post_partition_processing(
            input_file_path=text_input_path,
            best_partitions_by_cut_number_file=best_partitions_by_cut_number_file_path,
            cut_number=cut_number,
            positions_include_wt_aa_0idx=[mut['position'] for mut in params['mutations_0idx']],
            check_seq_complexity_idt=False,
            product_type="Gblock",
            validate_sample_number=10,
            validate_coding_start=None,
            min_dna_frag_length=params['provider_min_frag_len'],
            cost_per_nt=params['cost_per_nt']
        )
        log_dir = path.join(path.dirname(path.dirname(outfile_path)), 'logs')
        with open(path.join(log_dir, f'mutant_aa_fragments_{cut_number+1}fragments.json'), 'w', encoding='utf8') as f:
            json.dump(mutant_aa_fragments, f)
        with open(path.join(log_dir, f'mutant_dna_fragments_{cut_number+1}fragments.json'), 'w', encoding='utf8') as f:
            json.dump(mutant_dna_fragments, f)
        # copy the input used to generate the results into the result folder
        dst_path = path.join(path.dirname(outfile_path), f'{str(date.today())}_{path.basename(input_table_path)}')
        shutil.copyfile(src=input_table_path, dst=dst_path)
        outfile_paths.append(outfile_path)

    print(f'\n\033[1m'
          f'\n=================================================================================================='
          f'\n                        GENERATE & OPTIMIZE READY-TO-CLICK MODULES '
          f'\n                                     RUN COMPLETED!'
          f'\n=================================================================================================='
          f'\033[0m')

    return outfile_paths


def assemble_modules_and_generate_robot_instruction(input_table_path: str, ready_to_click_modules_path: str) -> str:
    """
    Assembles modules and generates a plate mapping sheet for robotic processing.

    This function reads desired variants from the input Excel file, generates possible variants
    based on module names or specified mutations, and validates the generated plate mapping sheet
    against the provided parameters such as plate format and plasmid IDs.

    Args:
        input_table_path (str): The file path to the input Excel file containing desired variants.
        ready_to_click_modules_path (str): The file path to the generated modules for mapping.

    Returns:
        None
    """
    desired_variants_input = pd.read_excel(input_table_path, sheet_name='input_desired_variants', header=0,
                                           index_col=None)
    out_file_path_, inputs_dict = transform_input_excel_sheet_to_text_input(input_table_path)
    desired_variants = desired_variants_input.to_dict(orient='records')
    plate_mapper_params = desired_variants_input.loc[0, ['plate_format', 'start_plasmid_id']]
    plate_format = int(plate_mapper_params['plate_format'])
    if str(plate_mapper_params['start_plasmid_id']) == 'nan':
        start_plasmid_id = None
    else:
        start_plasmid_id = str(plate_mapper_params['start_plasmid_id'])

    if desired_variants[0]['mutations'] == 'all':
        module_names = list(load_module_sheet(module_sheet_path=ready_to_click_modules_path)['Sequence Name'])
        desired_variant_muts_list = generate_all_possible_variants_from_modules(module_names)
        desired_variant_names = None
    else:
        desired_variant_muts_list = [variant['mutations'].split(',') for variant in desired_variants]
        desired_variant_names = [variant['names'] for variant in desired_variants]

    out_file_path = make_and_validate_plate_mapping_sheet(
        desired_variant_muts_list=desired_variant_muts_list,
        desired_variant_names=desired_variant_names,
        fragment_sheet_path=ready_to_click_modules_path,
        plate_format=plate_format,
        aa_seq=inputs_dict['SEQUENCE'],
        backbone_len=inputs_dict['BACKBONE_SIZE'],
        enzyme=inputs_dict['ENZYME'],
        five_prime_dna=inputs_dict['DNA_5_PRIME'],
        three_prime_dna=inputs_dict['DNA_3_PRIME'],
        start_plasmid_id=start_plasmid_id
    )
    return out_file_path

