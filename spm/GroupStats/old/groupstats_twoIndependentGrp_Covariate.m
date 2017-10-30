% This script performs 1 group or 2 group fMRI analysis
% It first loads configuration file containing groupstats parameters
% A template configuration file can be found at
% /home/fmri/fmrihome/SPM/spm8_scripts/GroupStats/groupstats_config.m.template
%
% This group-stats scripts are compatible with both Analyze and NIFTI formats
% To use either format, change the data type in groupStats_config.m.template
%
% Right now, this analysis script only supports one group
% analysis (one sample t test) and two group analysis (two sample t test).
%
% To run group analysis: type at Matlab command line:
% >> groupstats('groupstats_config.m')
%
% _________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: groupstats.m 2010-01-02 $
% -------------------------------------------------------------------------

function groupstats_Teresa(Config_File)

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
Config_File = strtrim(Config_File);
if(exist(Config_File,'file')==0)
  fprintf('Cannot find the configuration file ... \n');
  diary off; return;
end
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

% -------------------------------------------------------------------------
% Read in parameters
% -------------------------------------------------------------------------
stats_path       = strtrim(paralist.stats_path);
parent_folder    = strtrim(paralist.parent_folder);
subjlist_file    = strtrim(paralist.subjlist_file);
stats_folder     = strtrim(paralist.stats_folder);
output_folder    = strtrim(paralist.output_folder);
reg_file         = strtrim(paralist.reg_file);
template_path    = strtrim(paralist.template_path);
disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;
% -------------------------------------------------------------------------
% Check folder and file
% -------------------------------------------------------------------------

if ~exist(template_path,'dir')
  disp('Template folder does not exist!');
end

numgroup     = length(subjlist_file);

if(numgroup > 2)
  disp('Maximum number of groups should only be TWO');
  diary off; return;
end

if(exist(subjlist_file{1}, 'file') == 0)
  fprintf('File does not exist: %s\n', subjlist_file{1});
  diary off; return;
end

if(numgroup == 2)
  if(exist(subjlist_file{2}, 'file') == 0)
    fprintf('File does not exist: %s\n', subjlist_file{2});
    diary off; return;
  end
end

% -------------------------------------------------------------------------
% Check covariate of interest
% -------------------------------------------------------------------------
if ~isempty(reg_file)
  numreg = length(reg_file);
  for i = 1:numreg
    if ~exist(reg_file{i}, 'file')
      fprintf('Covariates file does not exist: %s\n', reg_file{i});
      diary off; return;
    end
  end
  
  reg_names = cell(numreg,1);
  for i = 1:numreg
    [a b c d e] = regexpi(reg_file{i}, '(\w+)\.\w+$');
    reg_names{i} = e{1}{1};
    reg_file{i} = fullfile(currentdir,reg_file{i});
  end
  reg_vec = 1;
  
else
  reg_vec = [];
end

% Initialize the batch system
spm_jobman('initcfg');
delete(get(0,'Children'));

% -------------------------------------------------------------------------
% Load subject list, constrast file and batchfile
% -------------------------------------------------------------------------
subjects1 = ReadList(subjlist_file{1});
batchfile   = 'batch_1group.mat';
numsubg1 = length(subjects1);

if(numgroup == 2)
  subjects2 = ReadList(subjlist_file{2});
  batchfile   = 'batch_2group.mat';
  numsubg2 = length(subjects2);
end

parent_folder = ReadList(parent_folder);

sub_stats1 = cell(numsubg1,1);
if numgroup == 2
  sub_stats2 = cell(numsubg2,1);
end

if isempty(parent_folder{1})
  for subcnt1 = 1:numsubg1
    secondpart = subjects1{subcnt1}(1:2);
    if str2double(secondpart) > 96
      pfolder = ['19' secondpart];
    else
      pfolder = ['20' secondpart];
    end
    sub_stats1{subcnt1} = fullfile(stats_path, pfolder, subjects1{subcnt1}, ...
      'fmri', 'stats_spm8', stats_folder);
  end
  if numgroup == 2
    for subcnt2 = 1:numsubg2
      secondpart = subjects2{subcnt2}(1:2);
      if str2double(secondpart) > 96
        pfolder = ['19' secondpart];
      else
        pfolder = ['20' secondpart];
      end
      sub_stats2{subcnt2} = fullfile(stats_path, pfolder, subjects2{subcnt2}, ...
        'fmri', 'stats_spm8', stats_folder);
    end
  end
else
  pfolder = parent_folder{1};
  for subcnt1 = 1:numsubg1
    sub_stats1{subcnt1} = fullfile(stats_path, pfolder, subjects1{subcnt1}, ...
      'fmri', 'stats_spm8', stats_folder);
  end
  if numgroup == 2
    for subcnt2 = 1:numsubg2
      sub_stats2{subcnt2} = fullfile(stats_path, pfolder, subjects2{subcnt2}, ...
        'fmri', 'stats_spm8', stats_folder);
    end
  end
end

contrastfile = fullfile(sub_stats1{1}, 'contrasts.mat');
load(contrastfile);

mkdir(output_folder);
cd (output_folder);

for i = 1:2:numTContrasts
  conNum = sprintf('%03d', i);
  conName = contrastNames{i};
  if numTContrasts == 1
    invName = ['-', contrastNames{i}];
  else
    invName = contrastNames{i+1};
  end
  dirName = fullfile(pwd, [conNum 'T_' conName]);
  mkdir(dirName);
  cd(dirName);
  load(fullfile(template_path,batchfile));
  % -----------------------------------------------------------------------
  % One Group Analysis
  % -----------------------------------------------------------------------
  if(numgroup == 1)
    if ~isempty(reg_vec)
      for j = 1:numreg
        regs = load(reg_file{j});
        matlabbatch{1}.spm.stats.factorial_design.cov(j).cname = reg_names{j};
        matlabbatch{1}.spm.stats.factorial_design.cov(j).c     = regs(:,1);
        matlabbatch{1}.spm.stats.factorial_design.cov(j).iCFT  = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(j).iCC   = 1;
        clear regs;
      end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = {};
    for j=1:numsubg1
      file = fullfile(sub_stats1{j}, ['con_0' conNum '.img,1']);
      matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{j} = file;
    end
    if isempty(reg_vec)
      matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = conName;
      matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = invName;
    else
      % x, first of contrast vector, x is 1 only if reg_vec indicates no
      % regressors of interest.
      x = ~ismember(1, reg_vec);
      matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [x  reg_vec];
      matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [x -reg_vec];
    end
  else
    % --------------------------------------------------------------------
    % Two Group Analysis
    % --------------------------------------------------------------------
    if ~isempty(reg_vec)
      for j = 1:numreg
        regs = load(reg_file{j});
        matlabbatch{1}.spm.stats.factorial_design.cov(j).cname = reg_names{j};
        matlabbatch{1}.spm.stats.factorial_design.cov(j).c     = regs(:,1);
        matlabbatch{1}.spm.stats.factorial_design.cov(j).iCFT  = 1;
        matlabbatch{1}.spm.stats.factorial_design.cov(j).iCC   = 1;
        clear regs;
      end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = {};
    for j=1:numsubg1
      file = fullfile(sub_stats1{j}, ['con_0' conNum '.img,1']);
      matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{j,1} = file;
    end
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = {};
    for j=1:numsubg2
      file = fullfile(sub_stats2{j}, ['con_0' conNum '.img,1']);
      matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{j,1} = file;
    end
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ['(g1-g2)(' conName ')'];
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ['(g2-g1)(' conName ')'];
  end
  
  % SPM.mat path for factorial desgin, SPM estimation and contrast.
  matlabbatch{1}.spm.stats.factorial_design.dir = {};
  matlabbatch{1}.spm.stats.factorial_design.dir{1} = pwd;
  matlabbatch{2}.spm.stats.fmri_est.spmmat{1}      = fullfile(pwd, 'SPM.mat');
  matlabbatch{3}.spm.stats.con.spmmat{1}           = fullfile(pwd, 'SPM.mat');
  
  save(batchfile,'matlabbatch');
  % Run group stats batch
  spm_jobman('run', ['./' batchfile]);
  cd ../
end

fprintf('Changing back to the directory: %s \n', currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('fMRI Group Stats finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
cd(currentdir);
diary off;
delete(get(0,'Children'));
clear all;
close all;
end