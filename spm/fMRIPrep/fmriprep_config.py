'''
Config file for fmri prep, a unified fmri preprocessing pipeline from the Poldrack lab (https://github.com/poldracklab/fmriprep)

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python fmriprep_config.py

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
working_dir = 'work_funcOnly10-11-17'

# Do you want to run freesurfer? 'no' = no; 'yes' = yes
freesurfer = 'no'

# Do you want to use another skullstripping template [choices: poldrack_oasis,scsnl; default scsnl]
skull_strip_template = 'scsnl'

# Do you want to run only one task? [specify name of task, e.g. resting_state_1 if yes; default is all]
task_t = 'all'

# Do you want to run aroma? 'no' = no; 'yes' = yes
aroma = 'no'

# Specify kernel width [4,6,8,10]
smoothing = 6

# Do you want to only preprocess your anatomicals? 'no' = process func and anat; 'yes' = process only anat 
# BEST PRACTICES: Run anat_only = 'yes' if you don't know how good your structural are. Then run anat_only = 'no' once you decide on the pipeline (below)
anat_only = 'yes'

# Do you want to only preprocess your functionals? 'no' = process func and anat; 'yes' = process only func
# BEST PRACTICES: Still experimental; only choose this option if you don't have any anatomical T1 images. If they are poor quality anatomicals, preprocess with swar pipeline (doesn't coregister to T1w; this option is --dismiss-t1w in other places)
func_only = 'no'

# Which pipeline do you want? swar (no T1w coregistration) or swgcar (T1w coregistration)
pipeline = 'swar'

######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'fmri_prep'

# Where to save output files from each job
job_output = project_dir+'/Jobs/fmriprep-%j.out'

# Where to save error files from each job
job_error = project_dir+'/Jobs/fmriprep-%j.err'

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
		'freesurfer':freesurfer,
 		'aroma':aroma,
		'pipeline':pipeline,
		'skull_strip_template':skull_strip_template,
		'task_t':task_t,
		'smoothing':smoothing,
		'anat_only':anat_only,
                'func_only':func_only,
		'working_dir':working_dir}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'fmriprep')
