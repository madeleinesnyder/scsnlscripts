'''
Config file for contrast generation

1. To run this you must use a custom python environment. Run the following commands:
	$ export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
	$ source activate

2. You can now run this script in bash. It will submit a parallel job for
   each of your subjects.
	$ python contrast_generator_config.py

'''

######################
## Script Variables ##
######################

#-----------------FILL OUT ALL THREE VARIABLES APPROPRIATLY----------

#Top level of project directory
projectdir = '/oak/stanford/groups/menon/projects/chamblin/test_project3'

#How many Conditions do you have in a SINGLE session (should have same number of
#conditions in each session)?
numcontrasts = 7


#How many sessions/runs will this be looking at? Sessons must be SAME size.
numsessions = 1


#If you want to compaire WITHIN sessions set variable to 1 else set
#variable to zero
comparewithin = 1


#If you are running ArtRepair and do not want to have movement correction
#make the movement correction to 0, 
#If you do want to factor in movement components set variable to 6
movementcorrection = 6

# Define your contrasts
# Make sure the even numbered contrasts are the opposite of the
# odd numbered contrasts. (i.e. all-rest; rest-all)


# SET THIS. Names of the contrasts:
contrastnames = ['CC','CC_neg','CE','CE_neg','EC','EC_neg','EE','EE_neg','cashout','cashout_neg','pump_lose','pump_lose_neg','trash','trash_neg','CC-CE','CE-CC','CC-EC','EC-CC','CC-EE','EE-CC','CE-EC','EC-CE','CE-EE','EE-CE','EC-EE','EE-EC']


# SET THIS. Contrasts defined in numbers, just based on your conditions. 
# So each vector should be as long as the number of conditions.
# NOTE: Each vector should sum to 0 unless contrasting with rest state
# If you named every other contrast the reverse of the previous one be sure to 
# do the same when making your contrast matrices

contrast = [[1, 0, 0, 0, 0, 0, 0], #[c1 c2 c3 ...] according to order in task design
	    [-1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0],
	    [0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0],
            [0, 0, -1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0]
            [0, 0, 0, 0, -1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0]
            [0, 0, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, -1],
            [1, -1, 0, 0, 0, 0, 0],
            [-1, 1, 0, 0, 0, 0, 0],            [1, 0, -1, 0, 0, 0, 0],
            [-1, 0, 1, 0, 0, 0, 0],
            [1, 0, 0, -1, 0, 0, 0],
            [-1, 0, 0, 1, 0, 0, 0],
            [0, 1, -1, 0, 0, 0, 0],
            [0, -1, 1, 0, 0, 0, 0],
            [0, 1, 0, -1, 0, 0, 0],
            [0, -1, 0, 1, 0, 0, 0],
            [0, 0, 1, -1, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 0]]


######################
### Slurm Variables ##
######################

# The name for your job
job_name = 'contrast_generator'

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
sys.path.append('/oak/stanford/groups/menon/scsnlscripts/spm/Chris/Utility')
from config import run_script

config_args = {'projectdir':projectdir,
	       'numcontrasts':numcontrasts,
	       'numsessions':numsessions,
	       'comparewithin':comparewithin,
	       'movementcorrection':movementcorrection,
	       'contrastnames':contrastnames,
	       'contrast':contrast}

slurm_args = {'job_name':job_name,
	      'job_output':job_output,
	      'job_error':job_error,
	      'job_time':job_time,
	      'job_qos':job_qos,
	      'job_nnodes':job_nnodes,
	      'job_mem':job_mem,
	      'job_partition':job_partition}

run_script(config_args,slurm_args,'contrast_generator')

