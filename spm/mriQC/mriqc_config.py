'''
Config file for mriqc, a machine-learning-based pipeline for removing artifacts and evaluating data quality from fmriprep from the Poldrack lab (https://github.com/poldracklab/mriqc)

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python mriqc_config.py

'''

######################
## Script Variables ##
######################

# Your top level project directory (full path)
project_dir = '/oak/stanford/groups/menon/fMRI_PREP_MAIDEN_VOYAGE/'

# Please specify the raw data directory (full path)
raw_dir = '/oak/stanford/groups/menon/rawdata/scsnl/'

# Please specify the .csv file or a cell array of subjects (full path)
subjectlist = '/oak/stanford/groups/menon/fMRI_PREP_MAIDEN_VOYAGE/data/subjectlist/scsnl_subjectlist.csv'

# Please specify the .txt file with your runs (full path)
runlist = '/oak/stanford/groups/menon/projects/mcsnyder/fMRI_PREP_MAIDEN_VOYAGE/data/subjectlist/runlist.txt'

# Please specify the name of the working directory (only the folder name)
working_dir = 'mriqc_work'


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'mriqc'

# Where to save output files from each job
job_output = project_dir+'/Jobs/mriqc-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/mriqc-%j.err'

# The amount of time to allot for each job
job_time = '16:00:00'

# Job priority
job_qos = 'normal'

# Number of nodes to request
job_nnodes = 1

# Amount of memory to allocate per CPU
job_mem = 20000 

# Partition to run your job on
# Options are:
#   - owners: access to all nodes
#   - menon: access to our lab nodes
#   - normal: the normal queue
job_partition = 'owners'

'''
Do not touch under this line
'''
import sys
sys.path.append('/oak/stanford/groups/menon/scsnlscripts/spm/Utility')
from config import run_script


config_args = {'raw_dir':raw_dir,
		'project_dir':project_dir,
		'subjectlist':subjectlist,
		'runlist':runlist,
		'working_dir':working_dir}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'mriqc')
