#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands
module load matlab
matlab -nodisplay -r "addpath('/home/groups/menon/oak/scsnlscripts/spm/Utility');matlab_dependencies;"$1";exit;"
