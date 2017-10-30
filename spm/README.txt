Scripts fully up to date (including PID directory fix) and tested:

TODO: Convert to jobarray submission rather than individual jobs

Preprocessing
-----------------
config: PreprocessfMRI/preprocessfmri_config.py
script: PreprocessfMRI/preprocessfmri.m

Individual Stats
----------------
config: IndividualStats/individualstats_config.py
script: IndividualStats/individualstats.m

Contrast Generation
---------------
config: ContrastGenerator/Contrast_Generator_config.m.template
script: ContrastGenerator/Contrast_Generator.m

TaskDesign
-------------
config: TaskDesign/taskdesign_generation_config.py
script: TaskDesign/taskdesign_generation.m


Scripts tested and working with Slurm but not with new Directory Structure:

Functional Connectivity
---------------
JN 5/26/17
config: FunctionalConnectivity/Seed_WholeBrain/fconnect_config.py
script: FunctionalConnectivity/Seed_WholeBrain/fconnect.m

gPPI
--------------
JN 5/30/17
config: EffectiveConnectivity/gPPI/gPPI_config.py
script: EffectiveConnectivity/gPPI/scsnl_gPPI.m

Contrast Change
--------------
JN 6/6/2017
config: ContrastChange/contrastchange_config.py
script: ContrastChange/constractchange.m

Movement Exclusion
---------------
CH 5/24/2017
config: Movement_Stats/movement_exclusion_config.py
script: Movement_Stats/movement_exclusion.m

Scripts converted for Slurm submit without PID structure but not tested:

EffectiveConnectivity/PPI/Effconn_phipi_config.py
EffectiveConnectivity/PPI/Effconn_ppi_config.py
EffectiveConnectivity/PPI/Effconn_ppi_volrep_config.py
DICOMtoSPGR/dicomtospgr_config.py
RSAnalysis/rsa_group_config.py
RSAnalysis/rsa_individual_config.py
RSAnalysis/rsa_wholebrain_config.py
ROIAnalysis/roi_signallevel_config.py
DTI/preprocess/dti_acpc_config.py
DTI/preprocess/dti_dcm2nii_config.py
DTI/preprocess/dti_preprocess_config.py
DTI/tractography/dti_betweenroitrack_config.py
DTI/tractography/dti_niftiroi2dtiroi_config.py
DTI/tractography/dti_trackstats_config.py
DTI/dti_wholebraintrack_config.py



