'''
Config file for fmri preprocessing

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
rawdir = '/scratch/PI/menon/projects/chamblin/rawdata'

#Path to the top level of the project directory
projectdir = '/scratch/PI/menon/projects/chamblin/test_project3'

#-Subject list
subjectlist = '/scratch/PI/menon/projects/chamblin/text_lists/subject_list1.txt'

#-Session list
sessionlist = '/scratch/PI/menon/projects/chamblin/text_lists/session_list1.txt'

#-Edat folder name
edatfolder = 'behavioral'

#-Global onset column (set the onset of real data collection)
#-current only 1 column is supported
#-If not specified, It will be read as the first REAL onset from OnsetColumns below
globalonsetcolumn = 'Fixation.OnsetTime'

#-Global onset shift: add a constant time to the first onset in 
#-'GlobalOnsetColumn' above (can be positive or negative)
globalonsetshift = 0


#-Onset column name
onsetcolumn = 'Equation.OnsetTime';

#-Offset column name (empty if not used, e.g. duration column is used)
offsetcolumn = ''

#-Duration column name (if specified, it will overide offset-onset)
durationcolumn = ''

#-Preset duration if no duration column AND no onset and offset columns (in seconds)
#-make sure that you set it to empty if you use Duration colum or onset and offset columns
presetduration = 9.5
# presetduration = 6

tasktypecolumn = ['Condition', 'Equation.ACC','Fixation']

#-Task names in the task_design.m (if empty, task names will be read from
#-the 'TaskTypeColumn' above, but be careful about the task order)
taskname = ['Task_acc', 'Control_acc', 'Task_inacc','Control_inacc']

#-Matched to the 'TaskTypeColumn' above
matchtaskname = [['Addition', 1,'*'], ['Control', 1,'*'], ['Addition', 0,'*'],['Control', 0,'*']]

#-Session name (make it empty to use the names in session list specified above)
# paralist.SessionName = {'fixed_block_1', 'fixed_block_2'};
sessionname = ''

#-Task design file name (*.m)
taskdesignfilename = 'taskdesign.m'

#-Rest exist or not
restexist = 1

######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'preprocessing_unflipped'

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
	       'edatfolder':edatfolder,
	       'globalonsetcolumn':globalonsetcolumn,
	       'globalonsetshift':globalonsetshift,
	       'onsetcolumn':onsetcolumn,
	       'offsetcolumn':offsetcolumn,
	       'durationcolumn':durationcolumn,
	       'presetduration':presetduration,
	       'tasktypecolumn':tasktypecolumn,
	       'taskname':taskname,
	       'matchtaskname':matchtaskname,
	       'sessionname':sessionname,
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


