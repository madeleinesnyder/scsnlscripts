%-remove global siganl from scrub with motion and interpolations

clear all; close all; clc;

subjects = ReadList('/mnt/mandarin1/Tianwen_Home_Temp/Projects/MovementCorrection/SCSNLData/Scripts/subjlist_scsnl_mvmnt.txt');

WholeBrainMask = '/mnt/mandarin1/Tianwen_Home_Temp/Projects/Mask/brain_mask.nii';

num_subj = length(subjects);

current_dir = pwd;

Mask = spm_read_vols(spm_vol(WholeBrainMask));
Mask(Mask ~= 0) = 1;
InBrainIndex = find(Mask ~= 0);

for isubj = 1:num_subj
  fprintf('---> %s\n', subjects{isubj});

  resting_folder = 'resting_state_1';
  
  output_subj_fmri_data_dir = fullfile('/mnt/mandarin1/Tianwen_Home_Temp/Projects/MovementCorrection/SCSNLData/Data/large_sample', ['20', subjects{isubj}(1:2)], subjects{isubj}, ...
    'fmri', resting_folder, 'smoothed_spm8_gfm');
  
  subj_fmri_data_file = fullfile(output_subj_fmri_data_dir, 'fswarI_motion_scrub_interp.nii');
  subj_fmri_data_V = spm_vol(subj_fmri_data_file);
  subj_fmri_data_D = spm_read_vols(subj_fmri_data_V);  
  
  img_dim = size(subj_fmri_data_D);
  
  subj_fmri_data_D = reshape(subj_fmri_data_D, prod(img_dim(1:3)), img_dim(4));
  subj_fmri_data_D = subj_fmri_data_D';
  
  gs = mean(subj_fmri_data_D(:, InBrainIndex), 2);
  
  nscans = length(gs);
  
  regressor_matrix = [gs, ones(nscans, 1)];
  
  subj_fmri_data = (eye(nscans) - regressor_matrix*pinv(regressor_matrix'*regressor_matrix)*regressor_matrix')*subj_fmri_data_D;
  
  subj_fmri_data = reshape(subj_fmri_data', img_dim);
  
  V = subj_fmri_data_V(1);
  V.dt(1) = 16;
  V.private.dat.dtype = 'FLOAT32-LE';
  V.private.dat.dim = img_dim(1:3);
  
  nscans = size(subj_fmri_data, 4);
  for i = 1:nscans
    if i < 10
      FlName = fullfile(output_subj_fmri_data_dir, ['gfswarI_motion_scrub_interp_00', num2str(i), '.nii']);
    else
      if i < 100
        FlName = fullfile(output_subj_fmri_data_dir, ['gfswarI_motion_scrub_interp_0', num2str(i), '.nii']);
      else
        FlName = fullfile(output_subj_fmri_data_dir, ['gfswarI_motion_scrub_interp_', num2str(i), '.nii']);
      end
    end    
    V.fname = FlName;
    V.private.dat.fname = FlName;    
    spm_write_vol(V, squeeze(subj_fmri_data(:,:,:,i)));
  end
  cd(output_subj_fmri_data_dir);
  
  unix('/usr/local/fsl/bin/fslmerge -t gfswarI_motion_scrub_interp.nii gfswarI_motion_scrub_interp_*.nii');
  unix('/bin/rm -rf gfswarI_motion_scrub_interp_*.nii');
  
  unix(sprintf('unpigz -fq %s', fullfile(output_subj_fmri_data_dir, 'gfswarI_motion_scrub_interp.nii.gz')));
  cd(current_dir);
end