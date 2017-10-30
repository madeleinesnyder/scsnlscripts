#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands
module load system
module load singularity/2.3.1
echo 'ENTERED'

singularity run --rm \ 
 -v $1:/oak/stanford/groups/menon/rawdata/scsnl \
 -v $2:/oak/stanford/groups/menon/fMRIPREP_MAIDEN_VOYAGE/data \
 -v /home/groups/menon/lab_shared/fmriprep:/home/groups/menon/lab_shared/software/miniconda2/lib/python2.7/site-packages:ro \
 poldracklab/fmriprep:latest /oak/stanford/groups/menon/rawdata/scsnl /oak/stanford/groups/menon/fMRIPREP_MAIDEN_VOYAGE/data \
 participant
 $3 \
 -t $4 \
 -w /oak/stanford/groups/menon/fMRIPREP_MAIDEN_VOYAGE/data/imaging
