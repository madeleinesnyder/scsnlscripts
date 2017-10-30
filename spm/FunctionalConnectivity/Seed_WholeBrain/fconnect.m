% Script for computing functional connectivity (resting state + seed-based)
%
% A template configuration file can be found at
% /home/fmri/fmrihome/SPM/spm8_scripts/FunctionalConnectivity/fconnect_config.m.template
%
% To run conversion: type at Matlab command line:
% >> fconnect(config_file)
% _________________________________________________________________________
% 2009-2011 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: fconnect.m Tianwen Chen 2011-08-17 v1$
% -------------------------------------------------------------------------

function fconnect (ConfigFile)

% Show the system information and write log files
warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('fMRI GroupAnalysis start at %d/%02d/%02d %02d:%02d:%02d\n',c);
disp('==================================================================');
fname = sprintf('groupanalysis-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

currentdir = pwd;

% -------------------------------------------------------------------------
% Check group analysis configuration and load it if it exists
% -------------------------------------------------------------------------
if ~exist(fullfile(ConfigFile), 'file')
    fprintf('Error: cannot find the configuration file... \n')
    return;
end

config = load(ConfigFile);
clear ConfigFile;

% -------------------------------------------------------------------------
% Read in parameters
% -------------------------------------------------------------------------
data_type = strtrim(config.data_type);
raw_dir = strtrim(config.datadir);
non_year_dir = strtrim(config.nonyear_dir);
fmri_type = strtrim(config.fmri_type);
project_dir = strtrim(config.project_dir);
stats_folder = strtrim(config.stats_dir);
processed_dir = fullfile(project_dir,'results',fmri_type,'participants');
subjectlist = strtrim(config.subject);
session_dir = strtrim(config.sessionlist);
preprocessed_dir = strtrim(config.preprocessed_dir);
pipeline_type = strtrim(config.pipeline);
TR_value = config.TR;
bandpass_include = config.bandpass;
fl_value = config.fl;
fh_value = config.fh;
roi_dir = strtrim(config.roi_dir);
roi_list = strtrim(config.roi_list);
numtrunc = config.numtrunc;


%-FC starting directory
run_FC_dir = pwd;
% Create local folder holding temporary data
tmpName = tempname;
temp_dir = fullfile(run_FC_dir, tmpName);
if exist(temp_dir,'dir')
  unix(sprintf('/bin/rm -rf %s', temp_dir));
end

disp('-------------- Contents of the Parameter List --------------------');
disp(config);
disp('------------------------------------------------------------------');
clear config;

%-Construct additional parameters
subjects = ReadList(subjectlist);
sublength = length(subjects);

%-Initialize path parameters
sessionlink_dir = cell(sublength,1);
imagedir = cell(sublength,1);
mvmntdir = cell(sublength,2);
subject_dir = cell(sublength,1);
pfolder = cell(sublength,1);

%-Update path parameters
if isempty(non_year_dir)
  for i = 1:sublength
    pfolder{i} = ['20', subjects{i}(1:2)];
  end
else
  for i = 1:sublength
    pfolder{i} = non_year_dir;
  end
end

%-update roi list
if ~isempty(roi_list)
  ROIName = ReadList(roi_list);
  NumROI = length(ROIName);
  roi_file = cell(NumROI, 1);
  for iROI = 1:NumROI
    ROIFile = spm_select('List', roi_dir, ['^', ROIName{iROI}]);
    if isempty(ROIFile)
      error('Folder contains no ROIs');
    end
    roi_file{iROI} = fullfile(roi_dir, ROIFile);
  end
end



for cnt = 1:sublength
    %sessionlink_dir{cnt} = fullfile(raw_dir, pfolder{cnt}, ....
     %                               subjects{cnt}, 'fmri', session_dir);

    imagedir{cnt} = fullfile(project_dir,'/data/imaging/participants', pfolder{cnt}, subjects{cnt}, ...
                            'fmri', session_dir, preprocessed_dir);

    mvmntdir{cnt,1} = fullfile(raw_dir, pfolder{cnt}, subjects{cnt}, ...
                               'fmri', session_dir, 'unnormalized');

    mvmntdir{cnt,2} = fullfile(project_dir,'/data/imaging/participants',pfolder{cnt}, subjects{cnt}, ...
                               'fmri', session_dir, preprocessed_dir);

    subject_dir{cnt} = fullfile(processed_dir, pfolder{cnt}, ...
                                subjects{cnt}, 'fmri',session_dir);
end

for FCi = 1:sublength
  disp('----------------------------------------------------------------');
  fprintf('Processing subject: %s \n', subject_dir{FCi});
  if exist(temp_dir, 'dir')
    unix(sprintf('/bin/rm -rf %s', temp_dir));
  end
  mkdir(temp_dir);
  cd(imagedir{FCi});
  fprintf('Copy files from: %s \n', pwd);
  fprintf('to: %s \n', temp_dir);
  if strcmp(data_type, 'nii')
    unix(sprintf('/bin/cp -af %s %s', [pipeline_type, '*.nii*'], temp_dir));
    if exist('unused', 'dir')
      unix(sprintf('/bin/cp -af %s %s', fullfile('unused', [pipeline_type, '*.nii*']), temp_dir));
    end
 else
    unix(sprintf('/bin/cp -af %s %s', [pipeline_type, '*.img*'], temp_dir));
    unix(sprintf('/bin/cp -af %s %s', [pipeline_type, '*.hdr*'], temp_dir));
    if exist('unused', 'dir')
      unix(sprintf('/bin/cp -af %s %s', fullfile('unused', [pipeline_type, '*.img*']), temp_dir));
      unix(sprintf('/bin/cp -af %s %s', fullfile('unused', [pipeline_type, '*.hdr*']), temp_dir));
    end
  end
  cd(temp_dir);
  unix('gunzip -fq *');
  newpipeline_type = pipeline_type;

  %-Bandpass filter data if it is set to 'ON'
  if bandpass_include == 1
    disp('Bandpass filtering data ......................................');
    bandpass_final_SPM(2, fl_value, fh_value, temp_dir, pipeline_type, data_type);
    disp('Done');
    %-Prefix update for filtered data
    newpipeline_type = ['filtered', pipeline_type];
  end

  %-Step 1 ----------------------------------------------------------------
  %-Extract ROI timeseries
  disp('Extracting ROI timeseries ......................................');
  [all_roi_ts, roi_name] = extract_ROI_timeseries(roi_file, temp_dir, 1, ...
                                      0, newpipeline_type, data_type);
  all_roi_ts = all_roi_ts';

  % Total number of ROIs
  numroi = length(roi_name);

  %-Step 2 ----------------------------------------------------------------
  %-Extract global signnals
  disp('Extract global signals .........................................');
  org_global_ts = ExtractGlobalSignal(data_type, newpipeline_type, temp_dir);

  %-Truncate ROI and global timeseries
  all_roi_ts = all_roi_ts(numtrunc(1)+1:end-numtrunc(2), :);
  global_ts = org_global_ts(numtrunc(1)+1:end-numtrunc(2));
  %-Run through multiple ROIs
  for roicnt = 1:numroi
    rts = all_roi_ts(:,roicnt);

    %-STEP 3 --------------------------------------------------------------
    %-Save covariates for each ROI
    disp('Making .txt file with timeseries and global signal ...........');
    unix(sprintf('gunzip -fq %s', fullfile(mvmntdir{FCi,1}, 'rp_I*')));
    unix(sprintf('gunzip -fq %s', fullfile(mvmntdir{FCi,2}, 'rp_I*')));
    rp2 = dir(fullfile(mvmntdir{FCi,2}, 'rp_I*.txt'));
    rp1 = dir(fullfile(mvmntdir{FCi,1}, 'rp_I*.txt'));
    if ~isempty(rp2)
      mvmnt = load(fullfile(mvmntdir{FCi,2}, rp2(1).name));
    elseif ~isempty(rp1)
      mvmnt = load(fullfile(mvmntdir{FCi,1}, rp1(1).name));
    else
      fprintf('Cannot find the movement file: %s \n', subjects{FCi});
      cd(current_dir);
      diary off; return;
    end



    %-Demeaned ROI timeseries and global signals
    rts = rts - mean(rts)*ones(size(rts, 1), 1);
    global_ts = global_ts - mean(global_ts)*ones(size(global_ts, 1), 1);
    NumVolsKept = length(global_ts);
    mvmnt = mvmnt(numtrunc(1)+1:numtrunc(1)+NumVolsKept, :);
    CovMtx = [global_ts mvmnt ones(NumVolsKept, 1)];

    %-Regress out global signals from ROI timeseries
    rts = (eye(NumVolsKept) - CovMtx*pinv(CovMtx'*CovMtx)*CovMtx')*rts;

    %-Covariates
    ts  = [rts global_ts mvmnt];

    reg_ts_dir = fullfile(subject_dir{FCi}, roi_name{roicnt}, 'timeseries');
    if ~exist(reg_ts_dir, 'dir')
      mkdir(reg_ts_dir);
    end

    reg_ts_all = fullfile(reg_ts_dir, 'roi_global_mvmnt.txt');
    save(reg_ts_all,'ts','-ascii','-tabs')

    %-STEP 4 --------------------------------------------------------------
    disp('Creating FC directories per subject ............................');
    final_FC_dir = fullfile(subject_dir{FCi}, roi_name{roicnt}, stats_folder);
    mkdir(final_FC_dir);
    cd   (final_FC_dir);
    if exist(fullfile(pwd,'SPM.mat'), 'file');
      unix('/bin/rm -rf *');
    end

    %-STEP 5 --------------------------------------------------------------
    disp('Creating task_design.mat for FC ................................');
    task_design_FC(reg_ts_all, 1)


    %-STEP 6 --------------------------------------------------------------
    disp('Running FC .....................................................');
    TR_val = TR_value;
    stats_fmri_fconnect_noscaling(temp_dir, data_type, numtrunc, newpipeline_type, TR_val);
    cd(run_FC_dir)
  end
  % Delete temporary folder  
  unix(sprintf('/bin/rm -rf %s', temp_dir));
  
end

cd(currentdir);

c     = fix(clock);
disp('==================================================================');
fprintf('Functional Connectivity finishes at %d/%02d/%02d %02d:%02d:%02d\n',c);
disp('==================================================================');

diary off;
clear all;
close all;

end
