% This script produce ROI statistics
% It helps you figure out if a particular ROI is in fact changing between 
% conditions in your experiment. Statistics are:
% *Percent Signal Change
% *t-score average
% *t-score percent voxels activated
% *beta average 
% 
% How to run it:
% start matlab with 'ml7spm8' and in matlab command prompt, type:
% roi_signallevel(config_file)
%
%__________________________________________________________________________
% 2009 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: roi_signallevel.m 2009-07-15 $
% -------------------------------------------------------------------------

function roi_signallevel(Config_File)

warning('off', 'MATLAB:FINITE:obsoleteFunction')
disp(['Current directory is: ',pwd]);
c     = fix(clock);
disp('==================================================================');
fprintf('ROI Signal Level Analysis start at %d/%02d/%02d %02d:%02d:%02d\n',c);
disp('==================================================================');
fname = sprintf('roi_signallevel-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

currentdir = pwd;

% -------------------------------------------------------------------------
% Check existence of the configuration file
% -------------------------------------------------------------------------

if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  return;
end

addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/marsbar'));

Config_File = strtrim(Config_File);
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

% Read in parameters
participant_folder = strtrim(paralist.participant_folder);
subjlist_file      = strtrim(paralist.subjlist_file);
stats_folder       = strtrim(paralist.stats_folder);
roi_folder         = strtrim(paralist.roi_folder);
tscore_threshold   = paralist.tscore_threshold;
roi_result_folder  = strtrim(paralist.roi_result_folder);

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

if ~exist(participant_folder, 'dir')
  fprintf('Folder does not exist: %s \n', participant_folder);
  diary off;
  return;
end

% Directory containing *_roi.mat files
if ~exist(roi_folder, 'dir')
  fprintf('Folder does not exist: %s \n', roi_folder);
  diary off;
  return;
end

% Read the list of subjects folder
subjects = ReadList(subjlist_file);
numsub = length(subjects);

% default, measure sc using entire event duration as coded in task_design
event_duration          = []; 

%--------------------------------------------------------------------------  
% run through all the subjects...
for ithsubject = 1:numsub
  disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<');
  fprintf('Processing subject: %s ...... \n', subjects{ithsubject});
  sub_stats_dir = fullfile(participant_folder,subjects{ithsubject},stats_folder); 
  if ~exist(sub_stats_dir, 'dir')
    fprintf('Folder does not exist: %s \n', sub_stats_dir);
    cd(currentdir);
    diary off; return;
  end

  % get percent signal change
  [signalchange{ithsubject}] = roi_signalchange_onesubject_scsnl_mrl(roi_folder, sub_stats_dir, event_duration); 
  
  % get tscore average and percent voxels activated in ROI
  [tscore_average{ithsubject}, tscore_percent_voxels{ithsubject}] = roi_tscore_onesubject(roi_folder,sub_stats_dir,tscore_threshold);
  
  % get beta average in ROI
  [beta_average{ithsubject}] = roi_beta_onesubject(roi_folder,sub_stats_dir);

end % subjects

% make a folder to hold roi statistics
mkdir(roi_result_folder);
cd(roi_result_folder);

% get summary data and stats for percent signal change, 
signal = signalchange; % change to generic name before saving
[signal_means, signal_stderr, signal_stats] = roi_stats_activation(signal, [], []); % get stats for all ROIs and events
save signalchange signal signal_means signal_stderr signal_stats

% get summary data and stats for tscore_average 
signal = tscore_average; % change to generic name before saving
[signal_means, signal_stderr, signal_stats] = roi_stats_activation(signal, [], []); % get stats for all ROIs and events
save tscore_average signal signal_means signal_stderr signal_stats

% get summary data and stats for tscore_percent_voxels  
signal = tscore_percent_voxels; % change to generic name before saving
[signal_means, signal_stderr, signal_stats] = roi_stats_activation(signal, [], []); % get stats for all ROIs and events
save tscore_percent_voxels signal signal_means signal_stderr signal_stats

% get summary data and stats for tscore_average 
signal = beta_average; % change to generic name before saving
[signal_means, signal_stderr, signal_stats] = roi_stats_activation(signal, [], []); % get stats for all ROIs and events
save beta_average signal signal_means signal_stderr signal_stats

PrintROIResults('signalchange');
PrintROIResults('tscore_average');
PrintROIResults('tscore_percent_voxels');
PrintROIResults('beta_average');

disp('-----------------------------------------------------------------');
fprintf('Changing back to the directory: %s \n', currentdir);
cd(currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('ROI Signal Level Analysis finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
diary off;
clear all;
close all;

end
