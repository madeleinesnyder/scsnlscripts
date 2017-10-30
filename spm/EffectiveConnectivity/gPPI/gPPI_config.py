'''
Config file for gPPI

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a job for each
   of your subjects in parallel
	$ python gPPI_config.py

'''


######################
## Script Variables ##
######################

# Please specify your top level project directory
project_dir = '/home/groups/menon/oak/projects/jnichola/testproject'

# Please specify the parent directory for the static data
# For YEAR data structure, use the first option
# For NONE YEAR data structure, use the second option
parent_dir = ''
#parent_dir = 'UCSFBoxer'

# Please specify the subject list .txt file
subjectlist = 'subjectlist.txt'

# Please specify the stats folder name
stats_dir = 'standard_addition1_4condition_1run_swaor'

# Please specify the prefix (pipeline) for your preprocessed images (pipeline)
pipeline = 'swar'

# Please specify a .txt file with the .nii file(s) for the ROI(s)
roifile_list = '/scratch/PI/menon/projects/jnichola/2017/script_testing/roifile_list.txt'

# Please specify a .txt file with the name(s) for the ROI(s)
roiname_list = '/scratch/PI/menon/projects/jnichola/2017/script_testing/roiname_list.txt'

# Please specify the task to include
# tasks_to_include = [ '1', 'task1', 'task2', 'task3'] -> must exist in all sessions
# tasks_to_include = [ '0', 'task1', 'task2', 'task3'] -> does not need to exist in all sessions
tasks_to_include = ['1', 'large_correct', 'small_correct', 'large_incorrect', 'small_incorrect']

# Confound names: leave these as default unless you have reason to change them
confound_names = ['R1','R2','R3','R4','R5','R6']


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'gPPI'

# Where to save output files from each job
# NOTE: you must make this directory!
job_output = project_dir+'/Jobs/repjob-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/repjob-%j.err'

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
job_partition = 'normal'

'''
Do not touch under this line
'''
import sys
sys.path.append('/home/groups/menon/oak/scsnlscripts/spm/Utility')
from config import run_script


config_args = {'project_dir':project_dir,
		'parent_dir':parent_dir,
		'subjectlist':subjectlist,
		'stats_dir':stats_dir,
		'pipeline':pipeline,
		'roifile_list':roifile_list,
		'roiname_list':roiname_list,
		'tasks_to_include':tasks_to_include,
		'confound_names':confound_names}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'scsnl_gPPI')

