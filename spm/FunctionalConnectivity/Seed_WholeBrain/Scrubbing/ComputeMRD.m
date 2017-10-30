%-spike regression criteria
clear all; close all; clc;

subjects = ReadList('/mnt/mandarin1/Tianwen_Home_Temp/Projects/MovementCorrection/SCSNLData/Scripts/subjlist_scsnl_mvmnt.txt');
output_file = 'summary_scsnl_46subjects.mat';

num_subj = length(subjects);

SubjectID = cell(num_subj, 1);
MvmntFile = cell(num_subj, 1);
RD = cell(num_subj, 1);
MeanRD = zeros(num_subj, 1);
Num_FD_Thresh = zeros(5, num_subj);
SpikeIndex = cell(num_subj, 1);

for isubj = 1:num_subj
  fprintf('---> %s\n', subjects{isubj});

  resting_folder = 'resting_state_1';
 
  mvmnt_file = fullfile('/fs/musk1', ['20', subjects{isubj}(1:2)], subjects{isubj}, ...
    'fmri', resting_folder, 'smoothed_spm8/rp_I.txt');
  
  load(mvmnt_file);
  mv_data = rp_I(9:end, :);
  
  nscans = size(mv_data, 1);
  
  %-FD
  fd = zeros(nscans, 1);  % Mean square displacement in two scans
  for i = 2:nscans
    fd(i) = (mv_data(i,1) - mv_data(i-1,1))^2 +...
      (mv_data(i,2) - mv_data(i-1,2))^2 +...
      (mv_data(i,3) - mv_data(i-1,3))^2 +...
      65*65*(mv_data(i,4) - mv_data(i-1,4))^2 +...
      65*65*(mv_data(i,5) - mv_data(i-1,5))^2 +...
      65*65*(mv_data(i,6) - mv_data(i-1,6))^2;
    fd(i) = sqrt(fd(i));
  end
  
  SubjectID{isubj} = subjects{isubj};
  MvmntFile{isubj} = mvmnt_file;
  RD{isubj} = fd;
  MeanRD(isubj) = mean(fd);
  
  for ithresh = 1:5
    Num_FD_Thresh(ithresh, isubj) = sum(fd > ithresh/10);
  end
  
  subj_spike_index = find(fd > 0.2);
  SpikeIndex{isubj} = subj_spike_index(:);
end

save(output_file, 'SubjectID', 'MvmntFile', 'RD', 'MeanRD', 'Num_FD_Thresh', 'SpikeIndex');