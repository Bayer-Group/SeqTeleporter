#!/usr/bin/env python
# coding: utf-8

# ## IMPORT REQUIRED RESOURCES
# --------------

# In[1]:


"""Imports"""
import re
import sys
import os
from os.path import join


SCRIPT_DIR = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(os.path.dirname(SCRIPT_DIR))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
print(f'Here is: {SCRIPT_DIR}')


from seqteleporter.partitioner.compute_best_partitions import compute_best_partitions, show_input_seq_info, prepare_compute_best_partitions_params
from seqteleporter.utils.load_input_params import validate_input_params


# In[2]:


def run(input_file_path):

    print('\033[1m==================================================================================================\033[0m')
    print('                                    \033[1m RUN STARTED! \033[0m ')
    print('\033[1m==================================================================================================\033[0m')

    validate_input_params(input_file_path)
    show_input_seq_info(input_file_path)
    params = prepare_compute_best_partitions_params(input_file_path)
    best_partitions_by_cut_number = compute_best_partitions(**params)


    print(f'\n\033[1m=================================================================================================='
          f'\n                                     RUN COMPLETED!                                                      '
          f'\n                    Total elapsed time: {best_partitions_by_cut_number["total_elapsed_time"]} seconds'
          f'\n==================================================================================================\033[0m')


# In[3]:


def list_input_files(input_dirs):
    sel_input_files = []
    for input_dir in input_dirs:
        for item in os.listdir(join(SCRIPT_DIR,  input_dir)):
            full_path = join(SCRIPT_DIR,  input_dir, item)
            if os.path.isfile(full_path) and re.search(f'[.]txt$',item):
                sel_input_files.append(full_path)
    return sel_input_files


# ## AVAILABLE HOSTS & ASSEMBLY CONDITIONS
# --------------

# ## RUN & SHOW RESULTS

# In[4]:


input_dirs = ['fix_cplx_diff_muts_inputs_T4_01h_25C']

input_files = list_input_files(input_dirs=input_dirs)[:2]

for input_file in input_files:
    run(input_file)

