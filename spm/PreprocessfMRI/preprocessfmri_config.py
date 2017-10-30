'''
Config file for fmri preprocessing

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

# Path to the directory with your raw data
raw_dir = '/oak/stanford/groups/menon/rawdata/scsnl'

# Your top level project directory
project_dir = '/oak/stanford/groups/menon/projects/tianwenc/2017_fmriprep/'

# Your subject list with full path
subjectlist = '/oak/stanford/groups/menon/projects/tianwenc/2017_fmriprep/data/subjectlist/oddball_subjects.csv'

# Your run list with full path
runlist = '/oak/stanford/groups/menon/projects/tianwenc/2017_fmriprep/data/subjectlist/oddball_runlist.txt'

# TR value for scans. Must be a float
tr_val = 2.0

# If you are running on files that have
# already been processed
inputimgprefix = ''

# The preprocessing pipeline to complete
# - "v" is the 1st version and "o" is the 2nd
#   version of the VolRepair pipeline
# - Choose from: 'swar',  'swavr', 'swaor',
#                'swgcar',  'swgcavr', 'swgcaor',
#		 'swfar', 'swfavr', 'swfaor',
#		 'swgcfar', 'swgcfavr', 'swgcfaor'
pipeline = 'swar'

# Name of the skullstripped T1w volume
# 'watershed_spgr' for brains stripped with mri_watershed
SPGRfilename = 'watershed_spgr'

# Whether you are running on 'task' or 'rest' fmri data
fmritype = 'task'


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'preprocessing'

# Where to save output files from each job
job_output = project_dir+'/Jobs/preprocessfmri-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/preprocessfmri-%j.err'

# The amount of time to allot for each job
job_time = '06:00:00'

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

config_args = {'raw_dir':raw_dir,
	       'subjectlist':subjectlist,
	       'runlist':runlist,
	       'tr_val':tr_val,
	       'inputimgprefix':inputimgprefix,
	       'pipeline':pipeline,
	       'SPGRfilename':SPGRfilename,
	       'project_dir':project_dir,
	       'fmritype':fmritype}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'preprocessfmri')

