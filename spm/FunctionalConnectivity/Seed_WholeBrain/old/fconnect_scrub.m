%-Functional connectivity using Power's scrubbing method
%
%-Tianwen Chen, 04/26/2012
% 2009-2012 SCSNL

function fconnect_scrub (Config_File)

warning('off', 'MATLAB:FINITE:obsoleteFunction')

Config_File = strtrim(Config_File);

if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off;
  return;
end

Config_File = Config_File(1:end-2);
eval(Config_File);

raw_data_dir = strtrim(paralist.raw_dir);
output_dir = strtrim(paralist.results_dir);
subjectlist = strtrim(paralist.subjectlist_file);
imagefilter = strtrim(paralist.pipeline_type);
session_folder = strtrim(paralist.session_dir);
preprocessed_folder = strtrim(paralist.preprocessed_dir);
bandpass_on = paralist.bandpass_include;
fl = paralist.fl_value;
fh = paralist.fh_value;
ROI_dir = strtrim(paralist.roi_dir);
ROI_list = strtrim(paralist.roi_list);
NUMTRUNC = paralist.numtrunc_value;

%==========================================================================
%-Store the current directory
current_dir = pwd;

data_type    = 'nii';
%-Create local folder holding temporary data (a 'temp' directory)
run_FC_dir = current_dir;
temp_dir = fullfile(run_FC_dir, 'temp');

%-Construct additional parameters
subjects = ReadList(subjectlist);
sublength = length(subjects);

%-Initialize path parameters
imagedir = cell(sublength,1);
mvmntdir = cell(sublength,2);
subject_dir = cell(sublength,1);

for cnt = 1:sublength
  YearID = ['20', subjects{cnt}(1:2)];
  imagedir{cnt} = fullfile(raw_data_dir, YearID, subjects{cnt}, ...
    'fmri', session_folder, preprocessed_folder);

  mvmntdir{cnt,1} = fullfile(raw_data_dir, YearID, subjects{cnt}, ...
    'fmri', session_folder, 'unnormalized');

  mvmntdir{cnt,2} = fullfile(raw_data_dir, YearID, subjects{cnt}, ...
    'fmri', session_folder, preprocessed_folder);

  subject_dir{cnt} = fullfile(output_dir, YearID, subjects{cnt}, 'fmri', 'resting_state');
end

if ~isempty(ROI_list)
  ROIName = ReadList(ROI_list);
  NumROI = length(ROIName);
  roi_file = cell(NumROI, 1);
  for iROI = 1:NumROI
    ROIFile = spm_select('List', ROI_dir, ['^', ROIName{iROI}]);
    if isempty(ROIFile)
      error('Folder contains no ROIs');
    end
    roi_file{iROI} = fullfile(ROI_dir, ROIFile);
  end
end

contrastNames = cell(2,1);
contrastVecs = cell(2,1);

contrastNames{1} = 'PosConn';
contrastNames{2} = 'NegConn';
contrastVecs{1} = 1;
contrastVecs{2} = 2;
numTContrasts = 2;

%-Copy raw data to your temporary folder and run filtering
for FCi = 1:sublength
  disp('----------------------------------------------------------------');
  fprintf('Processing subject: %s \n', subjects{FCi});
  if exist(temp_dir, 'dir')
    unix(sprintf('/bin/rm -rf %s', temp_dir));
  end
  mkdir(temp_dir);
  cd(imagedir{FCi});
  fprintf('Copy files from: %s \n', pwd);
  fprintf('to: %s \n', temp_dir);
  if strcmp(data_type, 'nii')
    if exist([imagefilter, 'I.nii.gz'], 'file')
      unix(sprintf('/bin/cp -af %s %s', [imagefilter, 'I.nii.gz'], temp_dir));
    else
      if exist([imagefilter, 'I.nii'], 'file')
        unix(sprintf('/bin/cp -af %s %s', [imagefilter, 'I.nii'], temp_dir));
      else
        disp('Raw image file does not exist');
        continue;
      end
    end
  end
  cd(temp_dir);
  unix('gunzip -fq *');

  %-Bandpass filter data if it is set to 'ON'
  if bandpass_on == 1
    disp('Bandpass filtering data ......................................');
    bandpass_final_SPM(2, fl, fh, temp_dir, imagefilter, data_type);
    disp('Done');
    %-Prefix update for filtered data
    newimagefilter = ['filtered', imagefilter];
  end
  %
  ImgFile = spm_select('FPList', temp_dir, ['^', newimagefilter, '.*\.nii']);
  ImgFile = deblank(cellstr(ImgFile));
  NumOrgImg = length(ImgFile);
  unix(sprintf('/bin/rm -rf %s', fullfile(temp_dir, [imagefilter, '*.nii'])));
  for i = 1:NUMTRUNC(1)
    unix(sprintf('/bin/rm -rf %s', ImgFile{i}));
  end
  for i = NumOrgImg-NUMTRUNC(2)+1:NumOrgImg
    unix(sprintf('/bin/rm -rf %s', ImgFile{i}));
  end

  ImgFile(1:NUMTRUNC(1)) = [];
  ImgFile(end-NUMTRUNC(2)+1:end) = [];

  %-Regress out global signal and movements
  org_global_ts = ExtractGlobalSignal(data_type, newimagefilter, temp_dir);
  global_ts = org_global_ts(1:end);
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
    return;
  end
  mvmnt = mvmnt(NUMTRUNC(1)+1:NUMTRUNC(1)+length(global_ts), :);

  %-Scrub criteria
  ScrubIndex = PowerScrub(mvmnt, char(ImgFile));

  for i = 1:length(ScrubIndex)
    unix(sprintf('/bin/rm -rf %s', ImgFile{ScrubIndex(i)}));
  end
  ImgFile(ScrubIndex) = [];

  global_ts(ScrubIndex) = [];
  mvmnt(ScrubIndex,:) = [];

  %-Demeaned ROI timeseries and global signals
  NumVolsKept = length(global_ts);
  global_ts = global_ts - mean(global_ts)*ones(size(global_ts, 1), 1);
  CovMtx = [global_ts mvmnt ones(NumVolsKept, 1)];

  %-Regress out global signals from ROI timeseries
  InvCovMtx = eye(NumVolsKept) - CovMtx*pinv(CovMtx'*CovMtx)*CovMtx';

  %-Extract ROI timeseries
  [all_roi_ts, roi_name] = extract_ROI_timeseries(roi_file, temp_dir, 1, ...
    0, newimagefilter, data_type);
  % Extract timeseries from all ROIs
  all_roi_ts = InvCovMtx*all_roi_ts';
  % Total number of ROIs
  numroi = length(roi_name);

  %-Run through multiple ROIs
  for roicnt = 1:numroi
    rts = all_roi_ts(:,roicnt);
    final_FC_dir = fullfile(subject_dir{FCi}, roi_name{roicnt}, 'stats_spm8');
    if ~exist(final_FC_dir, 'dir')
      mkdir(final_FC_dir);
    end
    cd(final_FC_dir);
    if exist(fullfile(final_FC_dir, 'SPM.mat'), 'file');
      unix('/bin/rm -rf *');
    end
    RegMtx = [rts, global_ts, mvmnt];
    RegFile = fullfile(final_FC_dir, 'roi_global_mvmnt.txt');
    save(RegFile, 'RegMtx', '-ascii', '-tabs')
    stats_fmri_fconnect_noscaling_noAR(temp_dir, data_type, 0, newimagefilter, RegFile);
    FlName = fullfile(final_FC_dir, 'contrasts.mat');
    save(FlName, 'contrastNames', 'contrastVecs', 'numTContrasts');
    cd(current_dir);
  end
end

cd(current_dir);

disp('Done');
end
