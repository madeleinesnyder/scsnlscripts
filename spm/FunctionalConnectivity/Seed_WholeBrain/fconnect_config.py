'''
Config file for functional connectivity

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a job for each
   of your subjects in parallel
	$ python fconnect_config.py

'''

######################
## Script Variables ##
######################

# Please specify the data type
# Default is 'nii'
data_type = 'nii'

# Path to the directory with your raw data
rawdir = '/home/groups/menon/oak/rawdata/scsnl'

# Please specify the data structure
# If it is a year structure:
nonyear_dir = ''
# If it is a non-year structure:
#nonyear_dir = 'foldercontainingsubjects'

# Please specify your top level project directory
project_dir = '/home/groups/menon/oak/projects/jnichola/testproject'

# Please specify the subject list .txt file
subjectlist = 'subjectlist.txt'

# Please specify the session name
session = 'resting_state'

# Please specify the preprocessed directory
preprocessed_dir = 'smoothed_spm8'

# Please specify the stats directory
stats_dir = 'stats_spm8'

# Please specify the prefix (pipeline) for your preprocessed images (pipeline)
pipeline = 'swar'

# Please specify the TR for your data (in seconds)
TR = 2.0

# Please specify whether to use bandpass filtering
# Set to 1 to bandpass filter, 0 to skip (default is 0)
bandpass = 1

# Please specify the bandpass filter parameters
# If you are not bandpassing then these values are ignored
# Lower frequency bound for filtering (in Hz)
# default is 0
fl = 0.008
# Upper frequency bound for filtering (in Hz)
# default is 0.2
fh = 0.1

# Please specify the folder containing your ROIs
roi_dir = '/scratch/PI/menon/projects/jnichola/2017/script_testing/rois'

# Please specify the ROI list
roi_list = '/scratch/PI/menon/projects/jnichola/2017/script_testing/roi_list.txt'

# Please specify the number of truncated images from the beginning and end
# (unit in SCANS not seconds, a two element vector, 1st slot for the beginning,
# and 2nd slot for the end, 0 means no truncation)
# Default is [0, 0]
numtrunc = [5, 5]

# Type of fmri you are running this on (either taskfmri or restfmri)
fmri_type = 'taskfmri'


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'fconnect'

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

config_args = {'data_type':data_type,
		'rawdir':rawdir,
		'nonyear_dir':nonyear_dir,
		'project_dir':project_dir,
		'subjectlist':subjectlist,
		'session':session,
		'preprocessed_dir':preprocessed_dir,
		'stats_dir':stats_dir,
		'pipeline':pipeline,
		'TR':TR,
		'bandpass':bandpass,
		'fl':fl,
		'fh':fh,
		'roi_dir':roi_dir,
		'roi_list':roi_list,
		'numtrunc':numtrunc,
		'fmri_type':fmri_type}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'fconnect')

