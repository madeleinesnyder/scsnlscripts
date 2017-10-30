#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands

module load system
module load singularity/2.3.1

#PROJECT_DIR=$1
#SUBJECT=$2
#VISIT=$3
#SESSION=$4

singularity exec -B /scratch:/scratch -B /oak/stanford/groups/menon:/oak/stanford/groups/menon/ -B /oak/stanford/groups/menon/software/singularity/dtiprep:/usr/local/miniconda/lib/python3.6/site-packages/dtiprep/ /share/PI/menon/lab_shared/singularity/singularity_images/X.img run_dwi_preproc $1 $2 $3 $4


