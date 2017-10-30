function scsnl_gPPI_volrep_mask(Config_File)
%-gPPI analysis for SCSNL data and analysis pipeline


%subject_list = '06-07-28.1';
%data_server = '/mnt/musk2';
%stats_folder = 'add_stats_spm8_swavr';
%roi_file_list = '/mnt/mabloo1/apricot1_share1/Longitudinal_TD_MD/addition_subtraction_td_md/connectivity/PPI/omnibus_peaks/right_ips/sphere_6-32_-60_44.nii';
%roi_name_list = 'right_ips';
% tasks_to_include = { '1', 'task1', 'task2', 'task3'} -> must exist in all sessions
% tasks_to_include = { '0', 'task1', 'task2', 'task3'} -> does not need to exist in all sessions
%tasks_to_include = {'1', 'complex', 'simple', 'find'};
%prep_pipeline = 'swavr';

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c = fix(clock);
disp('==================================================================');
fprintf('gPPI analysis started at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
fname = sprintf('scsnl_gPPI_volrep-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');
%==========================================================================
addpath(genpath('/home/fmri/fmrihome/SPM/spm8_scripts/Connectivity/EffectiveConnectivity/PPI/PPPI_v2012_1_22'));
addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/ArtRepair'));

%-Check existence of the configuration file
Config_File = strtrim(Config_File);

if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off;
  return;
end

Config_File = Config_File(1:end-2);

%-Read individual stats parameters
current_dir = pwd;
eval(Config_File);
clear Config_File;

%-Load parameters
data_server = strtrim(paralist.processed_dir);
subjects = ReadList(strtrim(paralist.subjectlist_file));
stats_folder = strtrim(paralist.stats_dir);
num_subj = length(subjects);
roi_file = ReadList(paralist.roi_file);
roi_name = ReadList(paralist.roiname_type);
num_roi_name = length(roi_name);
num_roi_file = length(roi_file);
tasks_to_include = paralist.tasks_list;
mask_file = paralist.mask_file;
confound_names = paralist.confounds_list;
prep_pipeline = strtrim(paralist.pipeline_type);

if num_roi_name ~= num_roi_file
  error('number of ROI files not equal to number of ROI names');
end

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

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

    year_id = ['20', subjects{i_subj}(1:2)];


    %-directory of SPM.mat file
    subject_stats_dir = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', stats_folder);

    subject_gPPI_stats_dir = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', [stats_folder, '_gPPI_mask']);

    subject_gPPI_stats_dir_temp = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', [stats_folder, '_gPPI_mask'], 'tempstats');

    if ~exist(subject_gPPI_stats_dir, 'dir')
      mkdir(subject_gPPI_stats_dir);
    else
      unix(sprintf('/bin/rm -rf %s', subject_gPPI_stats_dir));
      mkdir(subject_gPPI_stats_dir);
    end

    if ~exist(subject_gPPI_stats_dir_temp, 'dir')
      mkdir(subject_gPPI_stats_dir_temp);
    else
      unix(sprintf('/bin/rm -rf %s', subject_gPPI_stats_dir_temp));
      mkdir(subject_gPPI_stats_dir_temp);
    end

    cd(subject_gPPI_stats_dir_temp);

    unix(sprintf('/bin/cp -af %s %s', fullfile(subject_stats_dir, 'SPM.mat'), ...
      subject_gPPI_stats_dir_temp));

    unix(sprintf('/bin/cp -af %s %s', fullfile(subject_stats_dir, '*.img'), ...
      subject_gPPI_stats_dir_temp));

    unix(sprintf('/bin/cp -af %s %s', fullfile(subject_stats_dir, '*.hdr'), ...
      subject_gPPI_stats_dir_temp));


    P.subject = subjects{i_subj};
    P.directory = subject_gPPI_stats_dir_temp;

    %-Update the SPM path for gPPI analysis
    load('SPM.mat');
    SPM.swd = pwd;

    iG = [];
    col_name = SPM.xX.name;
    num_sess = numel(SPM.Sess);
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
      error('number of confound columns is larger than 6');
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
    img_name = cell(num_sess, 1);
    img_path = cell(num_sess, 1);
    num_scan = [1, SPM.nscan];

    for i_sess = 1:num_sess
      first_scan_sess = sum(num_scan(1:i_sess));
      img_name{i_sess} = SPM.xY.VY(first_scan_sess).fname;
      img_path{i_sess} = fileparts(img_name{i_sess});
      unix(sprintf('gunzip -fq %s', [img_name{i_sess}, '.gz']));
    end

    save SPM.mat SPM;
    clear SPM;

    %User input required (change analysis to be more specific)
    save(['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat'], 'P');
    PPPI_mask(['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat'], ['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat'], mask_file);
    gPPI_roi_dir = fullfile(subject_gPPI_stats_dir, ['PPI_', roi_name{i_roi}]);
    if ~exist(gPPI_roi_dir, 'dir')
      mkdir(gPPI_roi_dir);
    end
    scsnl_art_redo_mask(fullfile(subject_gPPI_stats_dir_temp, ['PPI_', roi_name{i_roi}]), prep_pipeline, gPPI_roi_dir, img_path, mask_file);

    for i_sess = 1:num_sess
      unix(sprintf('gzip -fq %s', img_name{i_sess}));
    end

    cd(subject_gPPI_stats_dir_temp);
    unix(sprintf('/bin/rm -rf %s', 'SPM.mat'));
    unix(sprintf('/bin/mv -f %s %s', '*.txt', gPPI_roi_dir));
    unix(sprintf('/bin/mv -f %s %s', '*.mat', gPPI_roi_dir));
    unix(sprintf('/bin/mv -f %s %s', '*.log', gPPI_roi_dir));

    unix(sprintf('/bin/rm -rf %s', subject_gPPI_stats_dir_temp));
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
