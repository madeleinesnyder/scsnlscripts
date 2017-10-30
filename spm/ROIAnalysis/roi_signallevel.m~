% This script produce ROI statistics
% It helps you figure out if a particular ROI is in fact changing between 
% conditions in your experiment. Statistics are:
% *Percent Signal Change
% *t-score average
% *t-score percent voxels activated
% *beta average 
%__________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% Tianwen Chen
% $Id: roi_signallevel.m rev.1 2010-01-24 $
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

%-Load configuration file
%-------------------------------------------------------------------------
if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  return;
end
Config_File = strtrim(Config_File);
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

%-Read in parameters
%--------------------------------------------------------------------------
project_dir        = strtrim(paralist.projectdir);
parent_dir         = strtrim(paralist.parent_dir);
subjectlist        = strtrim(paralist.subjectlist);
stats_dir          = strtrim(paralist.stats_dir);
roi_dir            = strtrim(paralist.roi_dir);
roilist            = strtrim(paralist.roilist);
tscore_threshold   = paralist.tscore_threshold;
roi_result_dir     = strtrim(paralist.roi_result_dir);
marsbar_path       = strtrim(paralist.marsbar_path);

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

%-Add marsbar to the search path
%--------------------------------------------------------------------------
if ~exist(marsbar_path, 'dir')
  fprintf('Marsbar toolbox does not exist: %s \n', marsbar_path);
  diary off;
  return;
else
  addpath(genpath(marsbar_path));
end

%-Check the roi_dir
%--------------------------------------------------------------------------
if ~exist(roi_dir, 'dir')
  fprintf('Folder does not exist: %s \n', roi_dir);
  diary off;
  return;
end

%-Read in subjects
%--------------------------------------------------------------------------
subjectlist = csvread(subjectlist,1);
subject = subjectlist(1);
subject = char(pad(string(subject),4,'left','0'));
numsub = length(subject);
visit = num2str(subjectlist(2));
session = num2str(subjectlist(3));

%-Construct the subject stats path
sub_stats = cell(numsub,1);
if isempty(parent_dir)
  for isubj = 1:numsub
    sub_stats{isubj} = fullfile(projectdir,subject,['visit',visit],['session',session],'fmri', 'stats_spm8', stats_dir);
  end
end

% default, measure sc using entire event duration as coded in task_design
event_duration          = []; 

%-ROI list
if ~isempty(roilist)
  ROIName = ReadList(roilist);
  NumROI = length(ROIName);
  roi_file = cell(NumROI, 1);
  for iROI = 1:NumROI
    ROIFile = spm_select('List', roi_dir, ROIName{iROI});
    if isempty(ROIFile) 
      error('Folder contains no ROIs'); 
    end
    roi_file{iROI} = fullfile(roi_dir, ROIFile);
  end
end
%--------------------------------------------------------------------------  
% run through all the subjects...
for ithsubject = 1:numsub
  disp('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<');
  fprintf('Processing subject: %s ...... \n', subject);
  sub_stats_dir = sub_stats{ithsubject};
  if ~exist(sub_stats_dir, 'dir')
    fprintf('Folder does not exist: %s \n', sub_stats_dir);
    cd(currentdir);
    diary off; return;
  end

  % get percent signal change
  [signalchange{ithsubject}] = roi_signalchange_onesubject_scsnl(roi_file, sub_stats_dir, event_duration); 
  
  % get tscore average and percent voxels activated in ROI
  [tscore_average{ithsubject}, tscore_percent_voxels{ithsubject}] = roi_tscore_onesubject(roi_file,sub_stats_dir,tscore_threshold);
  
  % get beta average in ROI
  [beta_average{ithsubject}] = roi_beta_onesubject(roi_file,sub_stats_dir);

end % subjects


% make a folder to hold roi statistics
mkdir(roi_result_dir);
cd(roi_result_dir);

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
