'''
Config file for fmri groupstats

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a job for groupstats
	$ python preprocessing_config.py

'''

######################
## Script Variables ##
######################

# Your top level project directory
project_dir = '/oak/stanford/groups/menon/projects/mcsnyder/PID_project'

# Type of fmri you are running this on (either taskfmri or restfmri)
fmri_type = 'taskfmri'

# Please specify the parent directory for the static data
# For YEAR data structure, use the first option
# For NONE YEAR data structure, use the second option
parent_dir = ''
#parent_dir = 'UCSFBoxer'

# Please specify the subject list .txt file
# For two group stats, specify two subject list files. e.g.: ['group1.txt','group2.txt']
subjectlist = '/oak/stanford/groups/menon/projects/mcsnyder/PID_project/data/subjectlist/subjectlist.csv'

# Please specify the stats folder name
stats_dir = 'stats_spm8'

# Please specify a folder name to put the analysis results
# Will be placed under project/results/groupstats
output_dir = 'groupstats'

# Please specify the file holding the regressors
# If there is no regressor, comment the first line and uncomment the second
#regfile = ['IQ.txt','EQ.txt']
regfile = ''

# Please specify the directory holding batch templates
template_dir   = '/oak/stanford/groups/menon/scsnlscripts/spm/BatchTemplates'


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'group_stats'

# Where to save output files from each job
job_output = project_dir+'/Jobs/groupstats-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/groupstats-%j.err'

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
import sys
sys.path.append('/oak/stanford/groups/menon/scsnlscripts/spm/Utility')
from config import run_script

config_args = {'project_dir':project_dir,
		'fmri_type':fmri_type,
		'parent_dir':parent_dir,
		'subjectlist':subjectlist,
		'stats_dir':stats_dir,
		'output_dir':output_dir,
		'regfile':regfile,
		'template_dir':template_dir}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'groupstats')

