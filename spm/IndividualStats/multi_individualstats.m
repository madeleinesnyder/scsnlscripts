function multi_individualstats (Config_File)

global currentdir idata_type run_imgdir template_path;

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('fMRI Multirun IndividualStats start at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
fname = sprintf('multi_individualstats-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

% -------------------------------------------------------------------------
% Check existence of the configuration file
% -------------------------------------------------------------------------

Config_File = strtrim(Config_File);

if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off;
  return;
end

% -------------------------------------------------------------------------
% Read individual stats parameters
% -------------------------------------------------------------------------
currentdir = pwd;
config = load(Config_File)
clear Config_File;

% Read in parameters
subject_i           = config.subject;
idata_type          = strtrim(config.image_type);
server_path         = strtrim(config.data_dir);
project_dir         = strtrim(config.project_dir);
fmri_type           = strtrim(config.fmri_type)
participant_path    = fullfile(project_dir,'results',fmri_type,'participants');
subjectlist         = strtrim(config.subjectlist);
exp_runlist         = {strtrim(config.runlist)};
task_dsgn           = strtrim(config.taskdesign_file);
pipeline            = strtrim(config.pipeline_type);
artpipeline         = strtrim(config.volpipeline_type);
contrastmat         = strtrim(config.contrast_file);
stats_folder        = strtrim(config.stats_dir);
template_path       = strtrim(config.template_dir);
preprocessed_folder = strtrim(config.preprocessed_dir);
include_mvmnt       = config.movement_include;
include_artrepair   = config.volrepair_include;
repaired_folder     = strtrim(config.volrepaired_dir);
repaired_stats      = strtrim(config.repairedstats_dir);

disp('-------------- Contents of the Parameter List --------------------');
disp(config);
disp('------------------------------------------------------------------');
clear config;

if ~exist(template_path,'dir')
  disp('Template folder does not exist!');
end

runs                      = ReadList(exp_runlist);
numruns                   = length(runs);
 
subjectlist_size          = size(subjectlist);
numgroup                  = subjectlist_size(1);
input_subjectlist_file    = subjectlist;
subjectlist = {};

for group_i = 1:numgroup
    subjectlist{group_i} = strtrim(input_subjectlist_file(group_i,:));
end

%subjectlist 		  = {subjectlist};
numpair 		  = length(subjectlist);
subjects 		  = cell(1,numpair);
sublength 		  = zeros(numpair,1);

for i = 1:numpair
  subjects{i} = subjectlist{i};
  sublength(i) = length(subjects{i});
end

if any(sublength ~= sublength(1))
  disp('Number of subjects are not the same across list files');
  cd(currentdir);
  diary off; return;
end
numsub = sublength(1);
clear sublength;

multi_pipeline = cell(numpair,1);
pipeline = ReadList(pipeline);

if length(pipeline) == 1
  for i = 1:numpair
    multi_pipeline{i} = pipeline{1};
  end
else
  multi_pipeline = pipeline;
end

multi_artpipeline = cell(numpair,1);
artpipeline = ReadList(artpipeline);

if length(artpipeline) == 1
  for i = 1:numpair
    multi_artpipeline{i} = artpipeline{1};
  end
else
  multi_artpipeline = artpipeline;
end

numrun = length(exp_runlist);
runs = cell(1,numpair);
if numrun > numpair
  disp('Number of runs is more than pairs of subjects');
  cd(currentdir);
  diary off; return;
end
if numrun > 1
  for i = 1:numrun
    temprun = ReadList(exp_runlist{i});
    if length(temprun) == 1
      runs{i} = repmat(temprun,numsub,1);
    else
      if length(temprun) ~= numsub
        disp('Number of experiments does not equal number of subjects');
        cd(currentdir);
        diary off; return;
      else
        runs{i} = temprun;
      end
    end
  end
else
  temprun = ReadList(exp_runlist);
  if length(temprun) == 1
    for i = 1:numpair
      runs{i} = repmat(temprun,numsub,1);
    end
  else
    for i = 1:numpair
      runs{i} = temprun;
    end
  end
end
clear temprun;

pfolder = cell(1,numpair);
run_dir = cell(1,numpair);
run_imgdir = cell(1,numpair);

for subcnt = 1:numsub
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
  fprintf('Processing subjects:');
  for pcnt = 1:numpair
    subjectlist   = csvread(input_subjectlist_file{pcnt},1);
    subject       = subjectlist(subject_i,1);
    subject       = char(pad(string(subject),4,'left','0'));
    visit         = num2str(subjectlist(subject_i,2));
    session       = num2str(subjectlist(subject_i,3));
    fprintf(' %s', subject);
  end
  fprintf('\n');
  stats_dir = fullfile(participant_path, subject,['visit',visit],['session',session], 'fmri', 'stats_spm8', stats_folder);
  if ~exist(stats_dir, 'dir')
    sprintf('Creating stats directory: %s \n', stats_dir);
    mkdir(stats_dir);
  else
    sprintf('Directory already exists: %s, deleting files inside \n', stats_dir);
    unix(sprintf('/bin/rm -rf %s', fullfile(stats_dir, '*')));
  end
  cd(stats_dir);
  for paircnt = 1:numpair
    run_dir{paircnt} = fullfile(project_dir, '/data/imaging/participants',subject,['visit',visit],['session',session], 'fmri', ...
                                 runs{paircnt}{subcnt});
    run_imgdir{paircnt} = fullfile(run_dir{paircnt}, preprocessed_folder);

    % Check the existence of preprocessed folder
    if ~exist(run_imgdir{paircnt}, 'dir')
      fprintf('Cannot find %s \n', run_imgdir{paircnt});
      cd(currentdir);
      diary off; return;
    end

    % If there is a ".m" at the end remove it.
    if(~isempty(regexp(task_dsgn, '\.m$', 'once' )))
      task_dsgn = task_dsgn(1:end-2);
    end
    % Load task_design file in raw server (SHOULD THIS BE IN DATA?)
    addpath(fullfile(run_dir{paircnt}));
    str = which(task_dsgn);
    if isempty(str)
      disp('Cannot find task design file in task_design folder.');
      cd(currentdir);
      diary off; return;
    end
    fprintf('Running the task design file: %s \n',str);
    eval(task_dsgn);
    rmpath(fullfile(run_dir{paircnt},'task_design'));
    % Unzip Image files in the preprocessed folder
    unix(sprintf('gunzip -fq %s', fullfile(run_imgdir{paircnt}, ...
                 [multi_pipeline{paircnt}, 'I*'])));

    % Update the design with the movement covariates
    if(include_mvmnt == 1)
      load task_design
      unix(sprintf('gunzip -fq %s', fullfile(run_imgdir{paircnt}, '*.txt.gz')));
      reg_file = spm_select('FPList', run_imgdir{paircnt}, '^rp_I');
      if isempty(reg_file)
        reg_path = fullfile(run_dir{paircnt}, 'unnormalized');
        unix(sprintf('gunzip -fq %s', fullfile(reg_path, '*.txt.gz')));
        reg_file = spm_select('FPList', reg_path, '^rp_I');
        if isempty(reg_file)
          disp('Cannot find the movement files');
          cd(currentdir);
          diary off; return;
        end
      end
      % Regressor names, ordered according regressor file structure
      reg_names = {'movement_x','movement_y','movement_z','movement_xr','movement_yr','movement_zr'};
      % 0 if regressor of no interest, 1 if regressor of interest
      reg_vec   = [0 0 0 0 0 0];
      disp('Updating the task design with movement covariates');
      save task_design.mat sess_name names onsets durations rest_exists reg_file reg_names reg_vec
    end
    if(numpair > 1)
      % Rename the task design file
      newtaskdesign = ['task_design_run' num2str(paircnt) '.mat'];
      movefile('task_design.mat', newtaskdesign);
    end
    % clear the variables used in input task_design.m file
    clear sess_name names onsets durations rest_exists reg_file reg_names reg_vec;
  end

  % Get the contrast file
  [pathstr, contrast_fname] = fileparts(contrastmat);
  if(isempty(pathstr) && ~isempty(contrast_fname))
    contrastmat = fullfile(currentdir,contrastmat);
  end
  if ~exist(contrastmat, 'file')
    disp('Cannot find contrast file');
    cd(currentdir);
    diary off; return;
  end

  % Run SPM multiple run individual stats
  individualfmri(multi_pipeline, numpair, contrastmat);

  % Redo analysis using ArtRepaired images and deweighting
  if include_artrepair == 1
    addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/ArtRepair'));
    repaired_folder_dir = cell(numpair, 1);
    for scnt = 1:numpair
      repaired_folder_dir{scnt} = fullfile(run_dir{scnt}, repaired_folder);
      unix(sprintf('gunzip -fq %s', fullfile(repaired_folder_dir{scnt}, ...
                     '*.txt.gz')));
      unix(sprintf('gunzip -fq %s', fullfile(repaired_folder_dir{scnt}, ...
                     [multi_artpipeline{scnt},'I*'])));
    end
    YearId = ['20', subjects{1}{subcnt}(1:2)];
    repaired_stats_dir = fullfile(participant_path, YearId, subjects{1}{subcnt}, 'fmri', 'stats_spm8', repaired_stats);
    if exist(repaired_stats_dir, 'dir')
      disp('------------------------------------------------------------');
      fprintf('%s already exists! Get deleted \n', repaired_stats_dir);
      disp('------------------------------------------------------------');
      unix(sprintf('/bin/rm -rf %s', repaired_stats_dir));
    end
    mkdir(repaired_stats_dir);
    scsnl_art_redo_multi(stats_dir, multi_artpipeline, repaired_stats_dir, ...
                   repaired_folder_dir);
    % copy contrasts.mat, task_design, batch_stats
    unix(sprintf('/bin/cp -af %s %s', fullfile(stats_dir, ['contrasts', '*']), ...
                 repaired_stats_dir));
    unix(sprintf('/bin/cp -af %s %s', fullfile(stats_dir, ['task_design', '*']), ...
                 repaired_stats_dir));
    unix(sprintf('/bin/cp -af %s %s', fullfile(stats_dir, 'batch_stats*'), ...
                 repaired_stats_dir));
    unix(sprintf('/bin/rm -rf %s', stats_dir));
    for scnt = 1:numpair
      unix(sprintf('gzip -fq %s', fullfile(repaired_folder_dir{scnt}, ...
                     [multi_artpipeline{scnt},'I*'])));
    end
  end
end

cd(currentdir);
% Change back to the directory from where you started.
fprintf('Changing back to the directory: %s \n', currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('fMRI Multirun Individual Stats finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
cd(currentdir);
diary off;
delete(get(0,'Children'));
clear all;
close all;

end

function individualfmri (multi_pipeline,numrun,contrastmat)

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
spm('defaults', 'fmri');
global idata_type run_imgdir template_path;

% Subject statistics folder
statsdir = pwd;

% -----------------------------------------------------------------------------
% fMRI design specification
% -----------------------------------------------------------------------------
load(fullfile(template_path,'batch_stats.mat'));

% Get TR value: initialized to 2 but will be update by calling GetTR.m
TR = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
% Initializing scans
matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = {};

for run = 1:numrun
  % Set preprocessed folder
  datadir = run_imgdir{run};

  %------------------------------------------------------------------------
  % Check the data type
  if isempty(idata_type)
    fselect = spm_select('List',datadir,['^',multi_pipeline{run},'I']);
    [strpath, fname, fext] = fileparts(fselect(1,:));
    if ismember(fext, {'.img', '.hdr'})
      data_type = 'img';
    else
      data_type = 'nii';
    end
  else
    data_type = idata_type;
  end
  %------------------------------------------------------------------------

  switch data_type
    case 'img'
      files = spm_select('ExtFPList', datadir, ['^',multi_pipeline{run},'I.*\.img']);
      nscans       = size(files,1);
    case 'nii'
      nifti_file = spm_select('ExtFPList', datadir, ['^',multi_pipeline{run},'I.*\.nii']);
      V       = spm_vol(deblank(nifti_file));
      nframes = V(1).private.dat.dim(4);
      files = spm_select('ExtFPList', datadir, ['^',multi_pipeline{run},'I.*\.nii'],1:nframes);
      nscans = size(files,1);
      clear nifti_file V nframes;
  end

  matlabbatch{1}.spm.stats.fmri_spec.sess(run) = ...
                                matlabbatch{1}.spm.stats.fmri_spec.sess(1);
  matlabbatch{1}.spm.stats.fmri_spec.sess(run).scans = {};

  % Input preprocessed images
  for nthfile = 1:nscans
    matlabbatch{1}.spm.stats.fmri_spec.sess(run).scans{nthfile} = ...
      deblank(files(nthfile,:));
  end


  if(numrun == 1)
    taskdesign_file = fullfile(statsdir, 'task_design.mat');
  else
    taskdesign_file = sprintf('%s/task_design_run%d.mat', statsdir, run);
  end

  reg_file = '';
  load(taskdesign_file);
  matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi{1}  = taskdesign_file;
  matlabbatch{1}.spm.stats.fmri_spec.sess(run).multi_reg = {reg_file};

end
matlabbatch{1}.spm.stats.fmri_spec.dir{1} = statsdir;

%--------------------------------------------------------------------------
% Estimation Setup
%--------------------------------------------------------------------------
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = strcat(statsdir,'/SPM.mat');

%--------------------------------------------------------------------------
% Contrast Setup
%--------------------------------------------------------------------------
matlabbatch{3}.spm.stats.con.spmmat{1} = strcat(statsdir,'/SPM.mat');

% Built the standard contrats only if the number of runs is one
% else use the user provided contrast file
if isempty(contrastmat)
  if (numrun >1 )
    disp(['The number of run is more than 1, No automatic contrast' ...
          ' generation option allowed, please spcify the contrast file']);
    diary off; return;
  else
    build_contrasts(matlabbatch{1}.spm.stats.fmri_spec.sess);
  end
else
  copyfile(contrastmat, './contrasts.mat');
end

load contrasts.mat;

for i=1:length(contrastNames)
  if (i <= numTContrasts)
    matlabbatch{3}.spm.stats.con.consess{i}.tcon.name   = contrastNames{i};
    matlabbatch{3}.spm.stats.con.consess{i}.tcon.convec = contrastVecs{i};
  elseif (i > numTContrasts)
    matlabbatch{3}.spm.stats.con.consess{i}.fcon.name = contrastNames{i};
    for j=1:length(contrastVecs{i}(:,1))
      matlabbatch{3}.spm.stats.con.consess{i}.fcon.convec{j} = ...
	  contrastVecs{i}(j,:);
    end
  end
end

save batch_stats matlabbatch
% Initialize the batch system
spm_jobman('initcfg');
delete(get(0,'Children'));
% Run analysis
spm_jobman('run', './batch_stats.mat');

for run = 1:numrun
  % Set scan data and stats directory
  datadir = run_imgdir{run};
  unix(sprintf('gzip -fq %s', fullfile(datadir, [multi_pipeline{run}, 'I*'])));
end

end
