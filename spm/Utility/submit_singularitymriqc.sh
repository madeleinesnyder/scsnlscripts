#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands

module load system
module load singularity/2.3.1

singularity exec -B /scratch:/scratch -B /oak/stanford/groups/menon:/oak/stanford/groups/menon/ /share/PI/menon/lab_shared/singularity/singularity_images/quality_check_poldrack.img mriqc $1 $2 participant -w $3  

