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

%-add warning message for two group stats with covariates, currently not
% supported yet to avoid confusion. - Tianwen, 12/11/2014
%
% To run group analysis: type at Matlab command line: 
% >> groupstats('groupstats_config.m')
%
% _________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: groupstats.m 2010-01-02 $
% -------------------------------------------------------------------------

function groupstats(ConfigFile)

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
stats_path       = strtrim(config.project_dir);
parent_dir       = strtrim(config.parent_dir);
subjlist_file    = strtrim(config.subjectlist);
stats_dir        = strtrim(config.statsfolder);
output_folder    = strtrim(config.outputfolder);
reg_file         = strtrim(config.regfile);
template_path    = strtrim(config.template_dir);
fmri_type        = strtrim(config.fmri_type);

disp('-------------- Contents of the Parameter List --------------------');
disp(config);
disp('------------------------------------------------------------------');
clear config;
% -------------------------------------------------------------------------
% Check folder and file
% -------------------------------------------------------------------------

output_folder = fullfile(stats_path,'/results/groupstats/',output_folder);
disp(output_folder)
stats_path = fullfile(stats_path,'/results/',fmri_type,'/participants/');

if ~exist(template_path,'dir')
  disp('Template folder does not exist!');
end

subjlist_size = size(subjlist_file);
numgroup = subjlist_size(1);
input_subjlist_file = subjlist_file;
subjlist_file = {};
for ii = 1:numgroup
    subjlist_file{ii} = strtrim(input_subjlist_file(ii,:));
end

regfile_size = size(reg_file);
numregs = regfile_size(1);
input_regfile = reg_file;
reg_file = {};
for ii = 1:numregs
   reg_file{ii} = strtrim(input_regfile(ii,:));
end

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
  if(numgroup == 2)
    fprintf('Covariate is not supported for two group stats\n');
    diary off; return;
  end
  numreg = length(reg_file);
  for i = 1:numreg
    if ~exist(reg_file{i}, 'file')
      fprintf('Covariates file does not exist: %s\n', reg_file{i});
      diary off; return;
    end
  end
  if(numgroup == 1)
    reg_names = cell(numreg,1);
    for i = 1:numreg
      [a b c d e] = regexpi(reg_file{i}, '(\w+)\.\w+$');
      reg_names{i} = e{1}{1};
      reg_file{i} = fullfile(currentdir,reg_file{i});
    end
    reg_vec = 1; 
  end
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
disp(subjects1)
batchfile   = 'batch_1group.mat';
numsubg1 = length(subjects1);
disp(numsubg1)

if(numgroup == 2)
  subjects2 = ReadList(subjlist_file{2});
  batchfile   = 'batch_2group.mat';
  numsubg2 = length(subjects2);
end

parent_dir = ReadList(parent_dir);

sub_stats1 = cell(numsubg1,1);
if numgroup == 2
  sub_stats2 = cell(numsubg2,1);
end
  
if isempty(parent_dir{1})
  for subcnt1 = 1:numsubg1
    secondpart = subjects1{subcnt1}(1:2);
    if str2double(secondpart) > 96
      pfolder = ['19' secondpart];
    else
      pfolder = ['20' secondpart];
    end
    sub_stats1{subcnt1} = fullfile(stats_path, pfolder, subjects1{subcnt1}, ...
                                   'fmri', 'stats_spm8', stats_dir);
  end
  if numgroup == 2
    for subcnt2 = 1:numsubg2
      sub_stats2{subcnt2} = fullfile(stats_path, pfolder, subjects2{subcnt2}, ...
                                     'fmri', 'stats_spm8', stats_dir);
    end
  end
  else
    pfolder = parent_dir{1};
    for subcnt1 = 1:numsubg1
      sub_stats1{subcnt1} = fullfile(stats_path, pfolder, subjects1{subcnt1}, ...
                                     'fmri', 'stats_spm8', stats_dir);
    end
    if numgroup == 2
      for subcnt2 = 1:numsubg2
        sub_stats2{subcnt2} = fullfile(stats_path, pfolder, subjects2{subcnt2}, ...
                                       'fmri', 'stats_spm8', stats_dir);
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
