#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands

# $1 is rawdata directory
# $2 is output directory
# $3 is trial_names
# $4 is onsetColname
# $5 is offsetColname
# $6 is durcolname
# $7 is presentdur
# $8 is tasktypecolumn
# $9 is subjectlist
# $10 is runlist

echo /oak/stanford/groups/menon/scsnlscripts/spm/TaskDesign/taskdesign_generation.py $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}

python /oak/stanford/groups/menon/scsnlscripts/spm/TaskDesign/taskdesign_generation.py $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
