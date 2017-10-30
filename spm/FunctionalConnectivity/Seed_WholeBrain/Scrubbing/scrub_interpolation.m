%-implement scrubbing method using temporal mask and interpolation
clear all; close all; clc;

subjects = ReadList('/mnt/mandarin1/Tianwen_Home_Temp/Projects/MovementCorrection/SCSNLData/Scripts/subjlist_scsnl_mvmnt.txt');

load('summary_scsnl_46subjects.mat');

num_subj = length(subjects);

current_dir = pwd;

for isubj = 1:num_subj
  fprintf('---> %s\n', subjects{isubj});

  resting_folder = 'resting_state_1';
 
  mvmnt_file = fullfile('/fs/musk1', ['20', subjects{isubj}(1:2)], subjects{isubj}, ...
    'fmri', resting_folder, 'smoothed_spm8/rp_I.txt');
  
  subj_fmri_data_dir = fullfile('/fs/musk1', ['20', subjects{isubj}(1:2)], subjects{isubj}, ...
    'fmri', resting_folder, 'smoothed_spm8');
  
  output_subj_fmri_data_dir = fullfile('/mnt/mandarin1/Tianwen_Home_Temp/Projects/MovementCorrection/SCSNLData/Data/large_sample', ['20', subjects{isubj}(1:2)], subjects{isubj}, ...
    'fmri', resting_folder, 'smoothed_spm8_gfm');
  
  exist_file = dir(fullfile(output_subj_fmri_data_dir, 'fswarI_motion_scrub_interp*'));
  if ~isempty(exist_file)
    unix(sprintf('/bin/rm -rf %s', fullfile(output_subj_fmri_data_dir, 'fswarI_motion_scrub_interp*')));
  end
  
  if ~exist(output_subj_fmri_data_dir, 'dir')
    mkdir(output_subj_fmri_data_dir);
  end
  
  unix(sprintf('unpigz -fq %s', fullfile(subj_fmri_data_dir, 'swarI.nii.gz')));
  
  subj_fmri_data_file = fullfile(subj_fmri_data_dir, 'swarI.nii');
  
  rp_I = load(mvmnt_file);
  
  %-linear detrend of movement parameters
  rp_I = spm_detrend(rp_I, 1);
  %-take out first 8 volumes
  mv_data = rp_I(9:end, :);
  
  nscans = size(mv_data, 1);
  
  deri_mvmnt = diff(rp_I);
  deri_mvmnt = deri_mvmnt(8:end, :);
  
  mvmnt_matrix = [mv_data, deri_mvmnt, mv_data.^2, deri_mvmnt.^2];
  %-construct movement matrix
  mvmnt_matrix = mvmnt_matrix - repmat(mean(mvmnt_matrix), nscans, 1);
  
  %-read swarI.nii
  subj_fmri_data_V = spm_vol(subj_fmri_data_file);
  subj_fmri_data_D = spm_read_vols(subj_fmri_data_V);  
  
  subj_fmri_data_D = subj_fmri_data_D(:,:,:,9:end);
  
  img_dim = size(subj_fmri_data_D);
  
  subj_fmri_data_D = reshape(subj_fmri_data_D, prod(img_dim(1:3)), img_dim(4));
  subj_fmri_data_D = subj_fmri_data_D';
  
  %-linear detrend of swarI.nii
  subj_fmri_data_D = spm_detrend(subj_fmri_data_D, 1);  
  
  %-frames with spikes
  subj_spike_index = SpikeIndex{isubj};
  num_spike = length(subj_spike_index);
  
  if num_spike ~= 0
    spike_regressor = zeros(nscans, num_spike);
    for ispike = 1:num_spike
      spike_regressor(subj_spike_index(ispike), ispike) = 1;
    end
    mvmnt_matrix(subj_spike_index, :) = 0;
  else
    spike_regressor = [];
  end
  
  regressor_matrix = [mvmnt_matrix, spike_regressor, ones(nscans, 1)];
  
  %-final spike regressors
  subj_fmri_data = (eye(nscans) - regressor_matrix*pinv(regressor_matrix'*regressor_matrix)*regressor_matrix')*subj_fmri_data_D;
  
  %-if there is spike, interpolate
  if num_spike ~= 0
    good_scans = 1:nscans;
    good_scans(subj_spike_index) = [];
    good_data = subj_fmri_data;
    good_data(subj_spike_index, :) = []; 
    interp_data = interp1(good_scans, good_data, subj_spike_index, 'nearest', 'extrap');
    subj_fmri_data(subj_spike_index, :) = interp_data;
    clear good_data interp_data;
  end
  
  subj_fmri_data = reshape(subj_fmri_data', img_dim);
  
  subj_fmri_data = BandPassProc(subj_fmri_data);
  
  if num_spike ~= 0 
    subj_fmri_data(:,:,:,subj_spike_index) = [];
  end
  
  V = subj_fmri_data_V(1);
  V.dt(1) = 16;
  V.private.dat.dtype = 'FLOAT32-LE';
  V.private.dat.dim = img_dim(1:3);
  
  nscans = size(subj_fmri_data, 4);
  for i = 1:nscans
    if i < 10
      FlName = fullfile(output_subj_fmri_data_dir, ['fswarI_motion_scrub_interp_00', num2str(i), '.nii']);
    else
      if i < 100
        FlName = fullfile(output_subj_fmri_data_dir, ['fswarI_motion_scrub_interp_0', num2str(i), '.nii']);
      else
        FlName = fullfile(output_subj_fmri_data_dir, ['fswarI_motion_scrub_interp_', num2str(i), '.nii']);
      end
    end    
    V.fname = FlName;
    V.private.dat.fname = FlName;    
    spm_write_vol(V, squeeze(subj_fmri_data(:,:,:,i)));
  end
  cd(output_subj_fmri_data_dir);
  
  unix('/usr/local/fsl/bin/fslmerge -t fswarI_motion_scrub_interp.nii fswarI_motion_scrub_interp_*.nii');
  unix('/bin/rm -rf fswarI_motion_scrub_interp_*.nii');
  
  unix(sprintf('pigz -fq %s', fullfile(subj_fmri_data_dir, 'swarI.nii')));
  cd(current_dir);
end



