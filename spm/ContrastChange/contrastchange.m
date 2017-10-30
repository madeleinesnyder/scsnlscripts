% This script performs contrast change for individual analysis results
% For new data structure
% _________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: contrastchange.m 2010-09-24 $
% -------------------------------------------------------------------------

function contrastchange (Config_File)

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('Individual Contrast Change starts at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
fname = sprintf('contrastchange-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

currentdir = pwd;

% -------------------------------------------------------------------------
% Check configuration and load it if it exists
% -------------------------------------------------------------------------
if ~exist(fullfile(Config_File), 'file')
    fprintf('Error: cannot find the configuration file... \n')
    return;
end

config = load(Config_File);
clear Config_File;

% Ignore white space if there is any
projectdir    = strtrim(config.projectdir);
participant_path =  fullfile(projectdir,'/results/taskfmri/participants/');
subjectlist         = strtrim(config.subject);
contrastmat         = strtrim(config.contrastfile);
stats_folder        = strtrim(config.statsdir);
template_path       = strtrim(config.template_dir);
disp('-------------- Contents of the Parameter List --------------------');
disp(config);
disp('------------------------------------------------------------------');
clear config;

if ~exist(participant_path, 'dir')
  disp('Cannot find the stats_server with individualstats ...');
  diary off; return;
end

if ~exist(contrastmat, 'file')
  fprintf('Cannot find contrast definition file ... \n');
  diary off; return;
end

subjects = ReadList(subjectlist);
numsub   = length(subjects);

load(contrastmat);
spm_jobman('initcfg');
delete(get(0,'Children'));

for subcnt = 1:numsub
  year_id = ['20', subjects{subcnt}(1:2)];
  subdir = fullfile(participant_path, year_id, subjects{subcnt}, ...
    'fmri', 'stats_spm8', stats_folder);
  cd(subdir);
  if exist(contrastmat,'file')
    delete(contrastmat);
  end
  if exist('batch_contrastchange.mat', 'file')
    delete('batch_contrastchange.mat');
  end
  condir = fullfile(currentdir,contrastmat);
  unix(sprintf('/bin/cp -af %s contrasts.mat', condir));
  load(fullfile(template_path, 'batch_contrastchange.mat'));
  matlabbatch{1}.spm.stats.con.spmmat = {};
  matlabbatch{1}.spm.stats.con.spmmat{1} = fullfile(subdir,'SPM.mat');
  matlabbatch{1}.spm.stats.con.delete = 1;
  for i=1:length(contrastNames)
    if (i <= numTContrasts)
      matlabbatch{1}.spm.stats.con.consess{i}.tcon.name   = contrastNames{i};
      matlabbatch{1}.spm.stats.con.consess{i}.tcon.convec = contrastVecs{i};
    elseif (i > numTContrasts)
      matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = contrastNames{i};
      for j=1:length(contrastVecs{i}(:,1))
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.convec{j} = ...
                                                      contrastVecs{i}(j,:);
      end
    end
  end
  save batch_contrastchange matlabbatch;
  clear matlabbatch;
  spm_jobman('run', './batch_contrastchange.mat');
end

fprintf('Changing back to the directory: %s \n', currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('Individual Contrast Change finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
cd(currentdir);
diary off;
delete(get(0,'Children'));
clear all;
close all;

end
