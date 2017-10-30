#!/bin/bash
#
#all SBATCH commands are given to this wrapper by the python config file
#below runs normal batch commands

module load system
module load singularity/2.3.1

echo "Singularity loaded, attempting to copy "$1 " to "$2

mkdir -p $2"/data/imaging/participants/"$4"/"

mkdir -p $2"/data/imaging/participants/"$4"/visit"$5"/"

mkdir -p $2"/data/imaging/participants/"$4"/visit"$5"/session"$6"/"

cp -r $1"/"$4"/visit"$5"/session"$6 $2"/data/imaging/participants/"$4"/visit"$5"/"
echo "SRC DIRECTORY"
echo $1"/"$4"/visit"$5"/session"$6
echo "TARGET DIRECTORY"
echo  $2"/data/imaging/participants/"$4"/visit"$5"/"

echo "Entering singularity shell with bindings (-B); executing fmriprep on participant "$3"; putting preprocessed data in "$2"/data/imaging/participants/"$4"/visit"$5"/session"$6"/task/preprocessed/"

singularity exec -B /scratch:/scratch -B /oak/stanford/groups/menon:/oak/stanford/groups/menon/ -B /oak/stanford/groups/menon/software/singularity/fmriprep_experimental/fmriprep:/usr/local/miniconda/lib/python3.6/site-packages/fmriprep/ /share/PI/menon/lab_shared/singularity/singularity_images/poldracklab_fmriprep_9-29-17.img fmriprep $2 $2 participant $3 -w $7 $8 $9 ${10} ${11} ${12} ${13} ${14} ${15}
