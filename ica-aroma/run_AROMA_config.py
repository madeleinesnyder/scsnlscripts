'''
#############################################################################
#################### Config file for run_AROMA.py ###########################
#############################################################################

To run:

1. To run this you must use a custom python 2 environment. It won't work in python 3! Run the following commands:
        $ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
        $ source activate scsnl

2. copy run_AROMA_config.py to your directory
3. fill out config and run: python run_AROMA_config.py

'''


######################
## Script Variables ##
######################

# Full path to your project directory
project_dir = '/oak/stanford/groups/menon/projects/jnichola/testproject'

# List of your task folders (should be a list not a txt file)
task_dirs = ['addition_block']

# 'nonaggr' to run the nonaggressive version of the algorithm and 'aggr' to run the aggressive version
aggr_sign = 'nonaggr'

# .csv file holding your subject list
subjectlist_file = '/oak/stanford/groups/menon/projects/jnichola/testproject/data/subjectlist/subjectlist.csv'

# specify which preliminary pipeline you have run (should be either 'swar' or 'swfar')
pipeline = 'swar'

# File with brain mask
brainmask = '/home/groups/menon/lab_shared/software/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'ica-aroma'

# Where to save output files from each job
job_output = project_dir+'/Jobs/ica-aroma-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/ica-aroma-%j.err'

# The amount of time to allot for each job
job_time = '01:00:00'

# Job priority
job_qos = 'normal'

# Number of nodes to request
job_nnodes = 1

# Amount of memory to allocate per CPU
job_mem = 4000

# Partition to run your job on
# Options are:
#   - owners: access to all nodes
#   - menon: access to our lab nodes
#   - normal: the normal queue
job_partition = 'owners'

'''
Do not touch under this line
'''

import sys, os
sys.path.append('/oak/stanford/groups/menon/scsnlscripts/ica-aroma/ICA-AROMA-master')
from ICA_AROMA_functions import submit_AROMA_job
import pandas as pd

slurm_args = {'job_name':job_name,
              'job_output':job_output,
              'job_error':job_error,
              'job_time':job_time,
              'job_qos':job_qos,
              'job_nnodes':job_nnodes,
              'job_mem':job_mem,
              'job_partition':job_partition}


subjectlist = pd.read_csv(subjectlist_file)
for i_row, row in subjectlist.iterrows():
    PID = str(row.PID).zfill(4)
    visit = str(row.visit)
    session = str(row.session)
    for task in task_dirs:
        subject_path = os.path.join(project_dir,'data/imaging/participants',PID,'visit'+visit,'session'+session,'fmri',task)
        submit_AROMA_job(subject_path,slurm_args,brainmask,pipeline,aggr_sign)








