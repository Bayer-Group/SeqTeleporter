#!/usr/bin/env python
# coding: utf-8

# ## IMPORT REQUIRED RESOURCES
# --------------

# In[2]:


"""Imports"""
import re
import sys
import os
from os.path import join


SCRIPT_DIR = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(os.path.dirname(SCRIPT_DIR))
sys.path.append(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
print(f'Here is: {SCRIPT_DIR}')


from seqteleporter.partitioner.compute_best_partitions import compute_best_partitions, get_lowest_cost_from_best_partitions, show_input_seq_info, prepare_compute_best_partitions_params
from seqteleporter.utils.utils import print_available_resources
from seqteleporter.utils.load_input_params import validate_input_params


# In[3]:


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


# In[4]:


def select_input(complexity,  input_dir):
    sel_input_files = []
    for item in os.listdir(join(SCRIPT_DIR,  input_dir)):
        full_path = join(SCRIPT_DIR,  input_dir, item)
        if os.path.isfile(full_path) and re.search(f'_{complexity}_',item):
            sel_input_files.append(full_path)
    return sel_input_files


# ## AVAILABLE HOSTS & ASSEMBLY CONDITIONS
# --------------

# In[5]:


print_available_resources()


# ## RUN & SHOW RESULTS

# In[9]:


input_dir = 'input_complexity_cost_rep0'
input_files = select_input(complexity=10000000000, input_dir=input_dir)

for input_file in input_files:
    input_file_path=join(SCRIPT_DIR, input_dir, input_file)
    run(input_file_path)

