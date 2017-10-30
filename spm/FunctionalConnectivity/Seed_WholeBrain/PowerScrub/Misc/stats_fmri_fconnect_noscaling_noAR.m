function stats_fmri_fconnect_noscaling_noAR (datadir, data_type, NUMTRUNC, imagefilter, RegMtx)

spm('Defaults', 'fmri');
statsdir = pwd;

% fMRI stats batchtemplate
load('/home/fmri/fmrihome/SPM/spm8_scripts/BatchTemplates/batch_stats.mat');

% Update TR
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2; 

% Select files (NIFTI/Analyze)
switch data_type
  case 'img'
    files = spm_select('ExtFPList', datadir, ['^', imagefilter,'.*\.img']);
    nscans       = size(files,1);            
  case 'nii'
    nifti_file = spm_select('ExtFPList', datadir, ['^', imagefilter,'.*\.nii']);
    V       = spm_vol(deblank(nifti_file));
    if length(V) == 4
      nframes = V.private.dat.dim(4);
      files = spm_select('ExtFPList', datadir, ['^', imagefilter,'.*\.nii'],1:nframes);
    else
      files = nifti_file;    
    end
    nscans = size(files, 1);
    clear nifti_file V nframes;
end

matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = {};
% Update input data
for filecnt = NUMTRUNC+1:nscans
  nthfile = filecnt - NUMTRUNC;
  matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{nthfile} = ...
    deblank(files(filecnt,:)); 
end

% Global scaling option and AR option
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';

% Update onsets
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi{1} = '';

% Update regressors
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {RegMtx};
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.dir{1} = statsdir;

% Update SPM.mat directory
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(statsdir, 'SPM.mat'); 

% Update contrasts
matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(statsdir,'SPM.mat');

matlabbatch{3}.spm.stats.con.consess{1}.tcon.name   = 'Positive Connectivity';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = 1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name   = 'Negative Connectivity';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = -1;

save batch_stats matlabbatch;

% Initialize the batch system
spm_jobman('initcfg');
delete(get(0,'Children'));
% Run analysis
spm_jobman('run', './batch_stats.mat');
delete(get(0,'Children'));
close all;

end