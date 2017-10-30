'''
Config file for fmri preprocessing

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python roi_signallevel_config.py

'''

######################
## Script Variables ##
######################
#
#__________________________________________________________________________
# 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
#
# $Id: roi_signallevel_config.m.template 2010-01-24 $
# -------------------------------------------------------------------------

# Please specify the server path
projectdir = '/oak/stanford/groups/menon/projects/yuanzh/Yuan_PNC'

# Please specify the parent folder
parent_dir = ''

# Please specify the subject list (full path of the .csv)
subjectlist  = '/home/groups/menon/oak/pncscripts/subjs_PNC_759.txt' # This should be a .csv for PID reading

# Please specify the folder containing SPM analysis results
stats_dir  = 'EI_D55'

# Please specify the folder (full path) holding defined ROIs
roi_dir    = '/home/groups/menon/oak/pncscripts/roi_folder/' # This should be changed to oak

# Please specify roi list file (full path)
roilist = '/home/groups/menon/oak/pncscripts/roifile_list.txt' # This should be changed to oak

# Please specify the t statistic threshold
tscore_threshold = 2.33

# Please specify the folder name to hold saved roi statistics
# You can change it to different studies or settings
roi_result_dir = 'roi_stats'

# Please specify the path of the Marsbar toolbox 
marsbar_path = '/home/groups/menon/lab_shared/software/spm8/toolbox/marsbar'


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'roi_signallevel'

# Where to save output files from each job
job_output = projectdir+'/Jobs/roi_signallevel-%j.out'

# Where to save error files from each job
job_error = projectdir+'/Jobs/roi_signallevel-%j.err'

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

config_args = {'projectdir':projectdir,
	       'subjectlist':subjectlist,
	       'parent_dir':parent_dir,
	       'subjectlist':subjectlist,
	       'stats_dir':stats_dir,
	       'roi_dir':roi_dir, 
               'roilist':roilist,
  	       'roi_result_dir':roi_result_dir,
	       'tscore_threshold':tscore_threshold,
	       'marsbar_path':marsbar_path}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'roi_signallevel')

