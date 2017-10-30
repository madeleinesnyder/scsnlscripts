'''
Config file for creating task design files 

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python taskdesign_config.py

'''

######################
## Script Variables ##
######################

# Path to the directory with your raw data (where edat files live)
raw_dir = '/oak/stanford/groups/menon/rawdata/scsnl/'

#Path to the top level of the project directory
project_dir = '/oak/stanford/groups/menon/projects/kduberg/SST_2017/'

#-Subject list
subjectlist = '/oak/stanford/groups/menon/projects/kduberg/SST_2017/data/subjectlist/test.csv'

#-Session list
runlist = '/oak/stanford/groups/menon/projects/kduberg/SST_2017/data/subjectlist/curr_run_list.txt'

#-Edat directory name
edat_dir = 'behavioral'

#-Global onset column (set the onset of real data collection)
#-current only 1 column is supported
#-If not specified, It will be read as the first REAL onset from OnsetColumns below
globalonsetcolumn = 'Starting.OnsetTime'

#-Global onset shift: add a constant time to the first onset in 
#-'GlobalOnsetColumn' above (can be positive or negative)
globalonsetshift = 0

#-Onset column name
onsetcolumn = 'GoSignalPhase.OnsetTime';

#-Offset column name ('none' if not used, e.g. duration column is used)
offsetcolumn = 'none'

#-Duration column name ('none' if not used, if specified, it will overide offset-onset)
durationcolumn = 'none'

#-Preset duration if no duration column AND no onset and offset columns (in seconds)
#-make sure that you set it to empty if you use Duration colum or onset and offset columns
presetduration = 9.5
# presetduration = 6

tasktypecolumn = 'stimulus_stop,GoSignalPhase.ACC,GoSignalPhase.RESP'

#-Task names in the task_design.m (if empty, task names will be read from
#-the 'TaskTypeColumn' above, but be careful about the task order)
taskname = 'go,trash,US,SS'

#-Matched to the 'TaskTypeColumn' above
matchtaskname = [[0, 1,'*'], [0, 0,'*'], [1, 1,''],[1, 0,'*'], [0, 0, '']]

#-Session name (make it empty to use the names in session list specified above)
# paralist.RunName = {'fixed_block_1', 'fixed_block_2'};
runname = ''

#-Task design file name (*.m)
taskdesignfilename = 'task_design.m'

#-Rest exist or not
restexist = 1

######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'task_design_generation'

# Where to save output files from each job
job_output = project_dir+'/Jobs/TaskDesign-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/TaskDesign-%j.err'

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
sys.path.append('/oak/stanford/groups/menon/scsnlscripts/spm/Utility')
from config import run_script

config_args = {'raw_dir':raw_dir,
	       'project_dir':project_dir,
	       'subjectlist':subjectlist,
	       'runlist':runlist,
	       'edat_dir':edat_dir,
	       'globalonsetcolumn':globalonsetcolumn,
	       'globalonsetshift':globalonsetshift,
	       'onsetcolumn':onsetcolumn,
	       'offsetcolumn':offsetcolumn,
	       'durationcolumn':durationcolumn,
	       'presetduration':presetduration,
	       'tasktypecolumn':tasktypecolumn,
	       'taskname':taskname,
	       'matchtaskname':matchtaskname,
	       'runname':runname,
	       'taskdesignfilename':taskdesignfilename,
	       'restexist':restexist}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'taskdesign_generation')


