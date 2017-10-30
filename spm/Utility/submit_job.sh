#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands
module load matlab

matlab -nodisplay -r "addpath('/oak/stanford/groups/menon/scsnlscripts/spm/Utility');matlab_dependencies; which $1;"$1";exit;"
