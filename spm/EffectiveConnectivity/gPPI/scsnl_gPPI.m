function scsnl_gPPI(ConfigFile)
%-gPPI analysis for SCSNL data and analysis pipeline

% Show the system information and write log files
warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('gPPI start at %d/%02d/%02d %02d:%02d:%02d\n',c);
disp('==================================================================');
fname = sprintf('gPPI-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

current_dir = pwd;

% -------------------------------------------------------------------------
% Check group analysis configuration and load it if it exists
% -------------------------------------------------------------------------
if ~exist(fullfile(ConfigFile), 'file')
    fprintf('Error: cannot find the configuration file... \n')
    return;
end

config = load(ConfigFile);
clear ConfigFile;

%-Load parameters
data_server = strtrim(config.project_dir);
data_server = fullfile(data_server,'/results/taskfmri/participants/');
stats_folder = strtrim(config.stats_dir);
parent_dir = strtrim(config.parent_dir);
subjects = {strtrim(config.subject)};
num_subj = 1;
roi_file = ReadList(config.roifile_list);
roi_name = ReadList(config.roiname_list);
num_roi_name = length(roi_name);
num_roi_file = length(roi_file);
tasks_to_include = config.tasks_to_include;
confound_names = config.confound_names;

% Python list --> Matlab cell str
tasks_size = size(tasks_to_include);
numtasks = tasks_size(1);
input_tasks = tasks_to_include;
tasks_to_include = {};
for ii = 1:numtasks
    tasks_to_include{ii} = strtrim(input_tasks(ii,:));
end

confound_size = size(confound_names);
numconfounds = confound_size(1);
input_cons = confound_names;
confound_names = {};
for ii = 1:numconfounds
    confound_names{ii} = strtrim(input_cons(ii,:));
end

if num_roi_name ~= num_roi_file
  error('number of ROI files not equal to number of ROI names');
end

for i_roi = 1:num_roi_file

  fprintf('===> gPPI for ROI: %s\n', roi_name{i_roi});

  load('ppi_master_template.mat');
  
  P.VOI = roi_file{i_roi};
  P.Region = roi_name{i_roi};
  P.extract = 'eig';
  P.CompContrasts = 0;
  P.Tasks = tasks_to_include;
  P.FLmask = 1;
  P.equalroi = 0;

  for i_subj = 1:num_subj

    fprintf('------> processing subject: %s\n', subjects{i_subj});

    if isempty(parent_dir)
        year_id = ['20', subjects{i_subj}(1:2)];
    else
        year_id = parent_dir;
    end

    %-directory of SPM.mat file
    subject_stats_dir = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', stats_folder);

    subject_gPPI_stats_dir = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', [stats_folder, '_gPPI']);

    if ~exist(subject_gPPI_stats_dir, 'dir')
      mkdir(subject_gPPI_stats_dir);
    end

    cd(subject_gPPI_stats_dir);

    unix(sprintf('/bin/cp -af %s %s', fullfile(subject_stats_dir, 'SPM.mat'), ...
      subject_gPPI_stats_dir));

    unix(sprintf('/bin/cp -af %s %s', fullfile(subject_stats_dir, '*.img'), ...
      subject_gPPI_stats_dir));

    unix(sprintf('/bin/cp -af %s %s', fullfile(subject_stats_dir, '*.hdr'), ...
      subject_gPPI_stats_dir));


    P.subject = subjects{i_subj};
    P.directory = subject_gPPI_stats_dir;

    %-Update the SPM path for gPPI analysis
    load('SPM.mat');
    SPM.swd = pwd;

    num_sess = numel(SPM.Sess);

    img_name = cell(num_sess, 1);
    img_path = cell(num_sess, 1);
    num_scan = [1, SPM.nscan];

    for i_sess = 1:num_sess
      first_scan_sess = sum(num_scan(1:i_sess));
      img_name{i_sess} = SPM.xY.VY(first_scan_sess).fname;
      img_path{i_sess} = fileparts(img_name{i_sess});
      unix(sprintf('gunzip -fq %s', [img_name{i_sess}, '.gz']));
    end


    iG = [];
    col_name = SPM.xX.name;

    num_confound = length(confound_names);

    for i_c = 1:num_confound
      iG_exp = ['^Sn\(.*\).', confound_names{i_c}, '$'];
      iG_match = regexpi(col_name, iG_exp);
      iG_match = ~cellfun(@isempty, iG_match);
      if sum(iG_match) == 0
        error('confound columns are not well defined or found');
      else
        iG = [iG find(iG_match == 1)];
      end
    end

    if length(iG) ~= num_confound*num_sess
      error('number of confound columns does not match SPM design');
    end

    num_col = size(SPM.xX.X, 2);
    FCon = ones(num_col, 1);
    FCon(iG) = 0;
    FCon(SPM.xX.iB) = 0;
    FCon = diag(FCon);

    num_con = length(SPM.xCon);

    %-make F contrast and run it
    SPM.xCon(end+1)= spm_FcUtil('Set', 'effects_of_interest', 'F', 'c', FCon', SPM.xX.xKXs);
    spm_contrasts(SPM, num_con+1);

    P.contrast = num_con + 1;

    SPM.xX.iG = sort(iG);
    for g = 1:length(iG)
      SPM.xX.iC(SPM.xX.iC==iG(g)) = [];
    end

    save SPM.mat SPM;
    clear SPM;
    %User input required (change analysis to be more specific)
    save(['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat'], 'P');
    
    PPPI(['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat']);

    for i_sess = 1:num_sess
      unix(sprintf('gzip -fq %s', img_name{i_sess}));
    end

    cd(subject_gPPI_stats_dir);
    unix(sprintf('/bin/rm -rf %s', 'SPM.mat'));
    unix(sprintf('/bin/rm -rf %s', '*.img'));
    unix(sprintf('/bin/rm -rf %s', '*.hdr'));
    unix(sprintf('/bin/mv -f %s %s', '*.txt', ['PPI_', roi_name{i_roi}]));
    unix(sprintf('/bin/mv -f %s %s', '*.mat', ['PPI_', roi_name{i_roi}]));
    unix(sprintf('/bin/mv -f %s %s', '*.log', ['PPI_', roi_name{i_roi}]));
  end
  cd(current_dir);
end

cd(current_dir);
disp('------------------------------------------------------------------');
fprintf('Changing back to the directory: %s \n', current_dir);
c     = fix(clock);
disp('==================================================================');
fprintf('gPPI analysis finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');

diary off;
clear all;
close all;
end
