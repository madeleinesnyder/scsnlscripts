'''
Config file for contrast change

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python contraschange_config.py

'''

######################
## Script Variables ##
######################

# Your top level project directory
projectdir = '/home/groups/menon/oak/projects/jnichola/testproject'

# Your subject list with full path
subjectlist = 'subjectlist.txt'

# Please specify the contrast definition file
contrastfile = 'new_contrast.mat'

# Please specify the stats folder that holds SPM.mat
statsdir = 'add_block'

# Please specify the directory holding batch templates (most likely leave this alone)
template_dir = '/home/groups/menon/oak/scsnlscripts/spm/BatchTemplates'

######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'contrastchange'

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
job_partition = 'owners'

'''
Do not touch under this line
'''
import sys
sys.path.append('/home/groups/menon/oak/scsnlscripts/spm/Utility')
from config import run_script

config_args = {'projectdir':projectdir,
	       'subjectlist':subjectlist,
	       'contrastfile':contrastfile,
	       'statsdir':statsdir,
	       'template_dir':template_dir}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'contrastchange')

