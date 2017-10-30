function scsnl_gPPI_volrep_empty_conds(Config_File)
%-gPPI analysis for SCSNL data and analysis pipeline

%stats_folder = 'add_stats_spm8_swavr';
%prep_pipeline = 'swavr';

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c = fix(clock);
disp('==================================================================');
fprintf('gPPI analysis started at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
fname = sprintf('scsnl_gPPI_volrep_empty_conds-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
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

if iscell(paralist.subjectlist_file)
    subjects = paralist.subjectlist_file;
else
    subjects = ReadList(strtrim(paralist.subjectlist_file));
end
stats_folder = strtrim(paralist.stats_dir);
num_subj = length(subjects);
roi_file = ReadList(paralist.roi_file);
roi_name = ReadList(paralist.roiname_file);
num_roi_name = length(roi_name);
num_roi_file = length(roi_file);
tasks_to_include = paralist.tasks_list;
confound_names = paralist.confounds_list;
prep_pipeline = strtrim(paralist.pipeline_type);

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

    year_id = ['20', subjects{i_subj}(1:2)];


    %-directory of SPM.mat file
    subject_stats_dir = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', stats_folder);

    subject_gPPI_stats_dir = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', [stats_folder, '_gPPI']);

    subject_gPPI_stats_dir_temp = fullfile(data_server, year_id, subjects{i_subj}, ...
      'fmri/stats_spm8', [stats_folder, '_gPPI'], 'tempstats');

    if ~exist(subject_gPPI_stats_dir, 'dir')
      mkdir(subject_gPPI_stats_dir);
    end

    if ~exist(subject_gPPI_stats_dir_temp, 'dir')
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
    checkiG = 1;

    for i_c = 1:num_confound
      cfd_i = confound_names{i_c};
      iG_exp = ['^Sn\(.\).', cfd_i, '$'];
      iG_match = regexpi(col_name, iG_exp);
      iG_match = ~cellfun(@isempty, iG_match);
      if sum(iG_match) == 0
        error('confound columns are not found');
      elseif sum(iG_match) > 1
	warning('confound columns are not unique -- this is OK if you have multisession data')
	checkiG = 0;
      else
        iG = [iG find(iG_match == 1)];
      end
    end

    if checkiG && length(iG) ~= 6
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


    %img_name = SPM.xY.VY(1).fname;
    %img_path = fileparts(img_name);
    %img_path = cellstr(img_path);

    save SPM.mat SPM;
    clear SPM;
    %unix(sprintf('gunzip -fq %s', [img_name, '.gz']));
    %User input required (change analysis to be more specific)
    save(['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat'], 'P');

    %JK CHANGES
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~HEADS UP ABOUT TO RUN SOME PPPI~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    addpath(genpath('/mnt/mapricot/musk2/home/fmri/fmrihome/SPM/spm8_scripts/Connectivity/EffectiveConnectivity/gPPI_empty_conds/'));
    which('PPPI')
    %JK CHANGES

    PPPI(['gPPI_', subjects{i_subj}, '_analysis_', roi_name{i_roi}, '.mat']);
    gPPI_roi_dir = fullfile(subject_gPPI_stats_dir, ['PPI_', roi_name{i_roi}]);
    if ~exist(gPPI_roi_dir, 'dir')
      mkdir(gPPI_roi_dir);
    end

    scsnl_art_redo(fullfile(subject_gPPI_stats_dir_temp, ['PPI_', roi_name{i_roi}]), prep_pipeline, gPPI_roi_dir, img_path);

%    unix(sprintf('gzip -fq %s', img_name));

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
%clear all;
close all;
end
