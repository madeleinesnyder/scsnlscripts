#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands
export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
source activate scsnl
echo $1
echo $2
echo $3
echo $4
python -s /oak/stanford/groups/menon/scsnlscripts/ica-aroma/run_AROMA.py $1 $2 $3 $4
