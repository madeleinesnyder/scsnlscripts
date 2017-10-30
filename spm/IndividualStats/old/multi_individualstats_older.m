function multi_individualstats (Config_File)

global currentdir idata_type sess_imgdir template_path;

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('fMRI Multisession IndividualStats start at %d/%02d/%02d %02d:%02d:%02d \n',c);
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

Config_File = Config_File(1:end-2);

% -------------------------------------------------------------------------
% Read individual stats parameters
% -------------------------------------------------------------------------
currentdir = pwd;
eval(Config_File);
clear Config_File;

% Read in parameters
idata_type          = strtrim(paralist.image_type);
server_path         = strtrim(paralist.raw_dir);
participant_path    = strtrim(paralist.processed_dir);
subjectlist         = strtrim(paralist.subjectlist_file);
exp_sesslist        = strtrim(paralist.sessionlist_file);
task_dsgn           = strtrim(paralist.taskdesign_file);
pipeline            = strtrim(paralist.pipeline_type);
artpipeline         = strtrim(paralist.volpipeline_type);
contrastmat         = strtrim(paralist.contrast_file);
stats_folder        = strtrim(paralist.stats_dir);
template_path       = strtrim(paralist.template_dir);
preprocessed_folder = strtrim(paralist.preprocessed_dir);
include_mvmnt       = paralist.movement_include;
include_artrepair   = paralist.volrepair_include;
repaired_folder     = strtrim(paralist.volrepaired_dir);
repaired_stats      = strtrim(paralist.repairedstats_dir);

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

if ~exist(template_path,'dir')
  disp('Template folder does not exist!');
end

numpair = length(subjectlist);
subjects = cell(1,numpair);
sublength = zeros(numpair,1);
for i = 1:numpair
  subjects{i} = ReadList(subjectlist{i});
  sublength(i) = length(subjects{i});
end;
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

numsess = length(exp_sesslist);
sessions = cell(1,numpair);
if numsess > numpair
  disp('Number of sessions is more than pairs of subjects');
  cd(currentdir);
  diary off; return;
end
if numsess > 1
  for i = 1:numsess
    tempsession = ReadList(exp_sesslist{i});
    if length(tempsession) == 1
      sessions{i} = repmat(tempsession,numsub,1);
    else
      if length(tempsession) ~= numsub
        disp('Number of experiments does not equal number of subjects');
        cd(currentdir);
        diary off; return;
      else
        sessions{i} = tempsession;
      end
    end
  end
else
  tempsession = ReadList(exp_sesslist);
  if length(tempsession) == 1
    for i = 1:numpair
      sessions{i} = repmat(tempsession,numsub,1);
    end
  else
    for i = 1:numpair
      sessions{i} = tempsession;
    end
  end
end
clear tempsession;

pfolder = cell(1,numpair);
sess_dir = cell(1,numpair);
sess_imgdir = cell(1,numpair);

for subcnt = 1:numsub
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
  fprintf('Processing subjects:');
  for pcnt = 1:numpair
    fprintf(' %s', subjects{pcnt}{subcnt});
  end
  fprintf('\n');
  YearId = ['20', subjects{1}{subcnt}(1:2)];
  stats_dir = fullfile(participant_path, YearId, subjects{1}{subcnt}, 'fmri', 'stats_spm8', stats_folder);
  if ~exist(stats_dir, 'dir')
    sprintf('Creating stats directory: %s \n', stats_dir);
    mkdir(stats_dir);
  else
    sprintf('Directory already exists: %s, deleting files inside \n', stats_dir);
    unix(sprintf('/bin/rm -rf %s', fullfile(stats_dir, '*')));
  end
  cd(stats_dir);
  for paircnt = 1:numpair
    pfolder{paircnt} = ['20',subjects{paircnt}{subcnt}(1:2)];
    sess_dir{paircnt} = fullfile(server_path, pfolder{paircnt}, ...
                                 subjects{paircnt}{subcnt}, 'fmri', ...
                                 sessions{paircnt}{subcnt});
    sess_imgdir{paircnt} = fullfile(sess_dir{paircnt}, preprocessed_folder);

    % Check the existence of preprocessed folder
    if ~exist(sess_imgdir{paircnt}, 'dir')
      fprintf('Cannot find %s \n', sess_imgdir{paircnt});
      cd(currentdir);
      diary off; return;
    end

    % If there is a ".m" at the end remove it.
    if(~isempty(regexp(task_dsgn, '\.m$', 'once' )))
      task_dsgn = task_dsgn(1:end-2);
    end
    % Load task_design file in raw server
    addpath(fullfile(sess_dir{paircnt}, 'task_design'));
    str = which(task_dsgn);
    if isempty(str)
      disp('Cannot find task design file in task_design folder.');
      cd(currentdir);
      diary off; return;
    end
    fprintf('Running the task design file: %s \n',str);
    eval(task_dsgn);
    rmpath(fullfile(sess_dir{paircnt}, 'task_design'));
    % Unzip Image files in the preprocessed folder
    unix(sprintf('gunzip -fq %s', fullfile(sess_imgdir{paircnt}, ...
                 [multi_pipeline{paircnt}, 'I*'])));

    % Update the design with the movement covariates
    if(include_mvmnt == 1)
      load task_design
      unix(sprintf('gunzip -fq %s', fullfile(sess_imgdir{paircnt}, '*.txt.gz')));
      reg_file = spm_select('FPList', sess_imgdir{paircnt}, '^rp_I');
      if isempty(reg_file)
        reg_path = fullfile(sess_dir{paircnt}, 'unnormalized');
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
    % Rename the task design file
    newtaskdesign = ['task_design_sess' num2str(paircnt) '.mat'];
    movefile('task_design.mat', newtaskdesign);
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

  % Run SPM multiple session individual stats
  individualfmri(multi_pipeline, numpair, contrastmat);

  % Redo analysis using ArtRepaired images and deweighting
  if include_artrepair == 1
    addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/ArtRepair'));
    repaired_folder_dir = cell(numpair, 1);
    for scnt = 1:numpair
      repaired_folder_dir{scnt} = fullfile(sess_dir{scnt}, repaired_folder);
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
fprintf('fMRI Multisession Individual Stats finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
cd(currentdir);
diary off;
delete(get(0,'Children'));
clear all;
close all;

end

function individualfmri (multi_pipeline,numsess,contrastmat)

% -------------------------------------------------------------------------
% Initialization
% -------------------------------------------------------------------------
spm('defaults', 'fmri');
global idata_type sess_imgdir template_path;

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

for sess = 1:numsess
  % Set preprocessed folder
  datadir = sess_imgdir{sess};

  %------------------------------------------------------------------------
  % Check the data type
  if isempty(idata_type)
    fselect = spm_select('List',datadir,['^',multi_pipeline{sess},'I']);
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
      files = spm_select('ExtFPList', datadir, ['^',multi_pipeline{sess},'I.*\.img']);
      nscans       = size(files,1);
    case 'nii'
      nifti_file = spm_select('ExtFPList', datadir, ['^',multi_pipeline{sess},'I.*\.nii']);
      V       = spm_vol(deblank(nifti_file));
      nframes = V(1).private.dat.dim(4);
      files = spm_select('ExtFPList', datadir, ['^',multi_pipeline{sess},'I.*\.nii'],1:nframes);
      nscans = size(files,1);
      clear nifti_file V nframes;
  end

  matlabbatch{1}.spm.stats.fmri_spec.sess(sess) = ...
                                matlabbatch{1}.spm.stats.fmri_spec.sess(1);
  matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans = {};

  % Input preprocessed images
  for nthfile = 1:nscans
    matlabbatch{1}.spm.stats.fmri_spec.sess(sess).scans{nthfile} = ...
      deblank(files(nthfile,:));
  end


  if(numsess == 1)
    taskdesign_file = fullfile(statsdir, 'task_design.mat');
  else
    taskdesign_file = sprintf('%s/task_design_sess%d.mat', statsdir, sess);
  end

  reg_file = '';
  load(taskdesign_file);
  matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi{1}  = taskdesign_file;
  matlabbatch{1}.spm.stats.fmri_spec.sess(sess).multi_reg = {reg_file};

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

% Built the standard contrats only if the number of sessions is one
% else use the user provided contrast file
if isempty(contrastmat)
  if (numsess >1 )
    disp(['The number of session is more than 1, No automatic contrast' ...
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

for sess = 1:numsess
  % Set scan data and stats directory
  datadir = sess_imgdir{sess};
  unix(sprintf('gzip -fq %s', fullfile(datadir, [multi_pipeline{sess}, 'I*'])));
end

end
