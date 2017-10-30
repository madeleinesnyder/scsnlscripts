#!/bin/bash
#
#SBATCH --job-name=prep_dwi
#SBATCH --time=06:30:00
#SBATCH --qos=normal
#SBATCH -p menon,owners,normal
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --mail-type=ALL

#ml biology fsl mrtrix

#export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
#source activate py27
#export PATH=${PATH}:/home/groups/menon/lab_shared/software/src/mrtrix3/bin
#export FREESURFER_HOME=/home/groups/menon/lab_shared/software/freesurfer6
#source /home/groups/menon/lab_shared/software/freesurfer6/SetUpFreeSurfer.sh
#export PATH=$PATH:/home/groups/menon/lab_shared/software/src/antsbin/bin
#ml load devel gcc

python /oak/stanford/groups/menon/scsnlscripts/spm/HARDI/dwi_preproc.py $1 $2 $3 $4

