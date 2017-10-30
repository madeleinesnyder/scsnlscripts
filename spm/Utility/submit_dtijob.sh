#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands

# $1 is rawdata directory
# $2 is output directory
# $3 is trial_names
# $4 is onsetColname

ml biology fsl mrtrix
export PATH=/home/groups/menon/lab_shared/software/miniconda2/bin:$PATH
source activate py27
export PATH=${PATH}:/home/groups/menon/lab_shared/software/src/mrtrix3/bin
export FREESURFER_HOME=/home/groups/menon/lab_shared/software/freesurfer6
source /home/groups/menon/lab_shared/software/freesurfer6/SetUpFreeSurfer.sh
export PATH=$PATH:/home/groups/menon/lab_shared/software/src/antsbin/bin
ml load devel gcc

/oak/stanford/groups/menon/scsnlscripts/spm/HARDI/run_dwi_preproc.sh $1 $2 $3 $4
