



Files
==========================================
MovementExclusion.py 
	python MovementExclusion.py (generates a .csv of the xyz roll pitch yaw stats across all subjects in the fmriprep directory)

plt_bold_confounds.py
	python ../../../scripts/prepfmri/plt_bold_confounds.py (plots the xyz roll pitch yaw and saves the figures for one subject. Run from fmriprep/sub-[subject]/func/ directory
	in progress 

fmriprep_config.py
	python fmriprep_config.py (see Pipieline Instruction section)



Pipeline Instruction
==========================================
python /oak/stanford/groups/menon/scsnlscripts/spm/fMRIPrep/fmriprep_config.py --> (run this with subjectlist and runlist)

squeue -u [sunetID] --> (use this command to check your job status)


Input/Output locations
==========================================
$OAK/projects/sunetID/[project_name]/Jobs --> contains all .err and .out files from slurm jobs submitted to the Sherlock nodes
.../data --> contains all raw data copied over from $OAK/rawdata/scsnl/, BIDS-formatted data (.../data/imaging/BIDS/), preprocessed folders
.../work --> contains all intermediate processing steps from the fmriprep pipeline, .html reports, good for debugging
.../fmriprep --> contains .html reports of preprocessing steps, the actual preprocessed data (anat and func), .tsv of func motion confounds
.../freesurfer --> contains fsaverage, fsaverage5, transforms, etc. 
.../scripts --> contains your specific fmriprep_config.py file, along with all other analysis .py and .m scripts copied from scsnlscripts/spm 


SLURM Specs
==========================================

partitions (-p)
---------------
dev 		Good for experimental code; gets jobs going quickly but gives the jobs fewer CPU and RAM resources
owners		Menon-owned nodes; plenty of CPU and RAM here, queueing dependent on how many lab users are on the nodes
normal		Open to all Sherlock users; slower resource allocation but may be faster than owners depending on useage
gpu 		For GPU processing (graphics-assisted aka deep-learning)
bigmem		Jobs requiring >64G RAM (probably not needed)

mem (--mem=X)
-------------
20G		fmriprep_config.py should use at least 20G for now

canceling/pausing jobs
---------------------
scancel -u sunetID	cancel all jobs
scancel <jobid>		cancel specific job
scontrol hold <jobid>	pause a job
scontrol resume <jobid> resume a job
control requeue <jobid> requeue a job (cancel and resume)


Under the Hood
==========================================
SCRIPTS
	fmriprep_config.py
	config.py
	submit_singularityfmriprep.sh
	
fmriprep_config.py passes [project_dir], [raw_dir], [subject_list.csv], [run_list.txt] and [SLURM arguments] to config.py
config.py constructs a string [callstr] and creates a custom submit_singularityfmriprep.sh script 
config.py call sbatch to execute the customized submit_singularityfmriprep.sh script with the slurm arguments for submitting the job
submit_singularityfmriprep.sh copies subject data from [raw_dir] to [project_dir] and generates a BIDS directory
submit_singularityfmriprep.sh opens a singularity container from the poldracklab docker image (converted from docker to singularity on scsnl PC)
submit_singularityfmriprep.sh executes the fmriprep command within that singularity container 
fmriprep preprocesses [project_dir] raw data and puts it in [project_dir] data/[subject_i]/[visit_i]/[session_i]/preprocessed/
fmriprep saves working results in [project_dir]/work


FAQ
==========================================
If the anatomical processing crashes will the functional processing still run?
	Yes. All the processes that don't depend on an anatomical will run.

If the anatomical processing crashes at a certain step, will the processing resume at the step if fmriprep is rerun?
	Yes but ONLY IF you specify the same working directory in the config file

Is there a flag to run only functional scans?
	No. But if you add functional scans to a directory where there is already a preprocessed anatomical image and a work directory with the anatomical processing steps completed, then the processed anatomical will be resued.


To debug:
******************************************

ENTER SINGULARITY SHELL WITH PROPER BINDINGS (-B) from docker image (.img)

	singularity shell \ 
	-B /lscratch:/scratch \ 
	-B /oak/stanford/groups/menon:/oak/stanford/groups/menon/ \
	-B /oak/stanford/groups/menon/software/singularity/fmriprep/fmriprep:/usr/local/miniconda/lib/python3.6/site-packages/fmriprep/ \
	/share/PI/menon/lab_shared/singularity/singularity_images/poldracklab_fmriprep_9-29-17.img 

WITHIN THE SHELL (you can't sudo...) 

	fmriprep [raw_dir] [project_dir] participant [subject_i] [visit_i] [session_i] -w [project_dir]/work  

EXAMINE ERROR LOG IN

	[project_dir]/Jobs
	




