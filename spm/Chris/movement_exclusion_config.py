'''
Config file for fmri movement exclusion

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python preprocessfmri_config.py

'''

######################
## Script Variables ##
######################

#raw directory: directory holding unnormalized data
rawdir = '/scratch/PI/menon/projects/chamblin/rawdata'

#-Top level of project directory
projectdir = '/scratch/PI/menon/projects/chamblin/test_project'

#-Folder name that contains rp_I.txt (by default, it is the preprocessed
# data folder
preprocessedfolder = 'smoothed_spm8'

#-Subject list file
subjectlist = '/scratch/PI/menon/projects/chamblin/text_lists/subject_list_unflipped.txt'

#-Session list file
sessionlist = '/scratch/PI/menon/projects/chamblin/text_lists/session_list1.txt'

#-Scan-to-scan threshold (unit in voxel)
scantoscancrit = 0.5;

######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'movement_exclusion'

# Where to save output files from each job
job_output = projectdir+'/Jobs/repjob-%j.out'

# Where to save error files from each job
job_error = projectdir+'/Jobs/repjob-%j.err'

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
job_partition = 'menon'

'''
Do not touch under this line
'''
import sys
sys.path.append('/home/groups/menon/oak/scsnlscripts/spm/Utility')
sys.path.append('/home/groups/menon/oak/scsnlscripts/spm/Chris')
from config import run_script

config_args = {'rawdir':rawdir,
	       'projectdir':projectdir,
	       'subjectlist':subjectlist,
	       'sessionlist':sessionlist,
	       'preprocessedfolder': preprocessedfolder,
	       'scantoscancrit':scantoscancrit}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'movement_exclusion')



