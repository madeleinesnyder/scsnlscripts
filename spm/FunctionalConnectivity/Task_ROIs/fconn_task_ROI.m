function fconn_task_ROI(Config_File)

% This script calculates the functional connectivity of ROI pairs
% Check the configurations in fconn_task_ROI_config.m.template
% Use new data structure
%
% To run, type 'ml7spm8' in unix, run fconn_task_ROI('config_file')
%__________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
% Tianwen Chen
% 03/36/10, v1
% 05/01/10, v2: writing results to txt files
% 09/15/10, v3: add option of removing global time series
%--------------------------------------------------------------------------

scsnl_id = '$Id: fconn_task_ROI.m 2010-09-15 v3 $';

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('ROI Functional Connectivity started at %d/%02d/%02d %02d:%02d:%02d \n',c);
fprintf('%s \n', scsnl_id);
disp('==================================================================');
fname = sprintf('fconn_task_ROI-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

%-Check existence of the configuration file
Config_File = strtrim(Config_File);

if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off; 
  return;
end

Config_File = Config_File(1:end-2);

%-Read individual stats parameters
currentdir = pwd;
eval(Config_File);
clear Config_File;

server_path   = strtrim(paralist.server_path);
parent_folder = strtrim(paralist.parent_folder);
subjectlist   = strtrim(paralist.subjectlist);
stats_folder  = strtrim(paralist.stats_folder);
roi_folder    = strtrim(paralist.roi_folder);
roi_list      = strtrim(paralist.roi_list);
global_ts     = paralist.rm_global_ts;
num_trunc     = paralist.num_trunc;
filt          = logical(paralist.filter);
TR            = paralist.TR;

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

%-Check existence of files
if ~exist(roi_folder, 'dir')
  fprintf('The ROI folder does not exist: %s \n', roi_folder);
  diary off; return;
end

%-Read in list files or cell array
parent_folder = multi_readlist(parent_folder);
numparent = length(parent_folder{1});
if numparent > 2
  disp('No more than two parent folder allowed');
  diary off; return;
end

subjects = multi_readlist(subjectlist);
numgroup = length(subjects);
numsub = zeros(numgroup,1);
for i = 1:numgroup
  numsub(i) = length(subjects{i});
end
if numgroup > 2
  disp('More than 2 groups!');
  diary off; return;
end

if numgroup == 2 && numparent ==1
  parent_folder{1}{2} = parent_folder{1}{1};
end

%-Start to process ROI functional connectivity
results_fconn_task_ROI = {};
for groupcnt = 1:numgroup
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
  fprintf('Processing Group %d \n', groupcnt);
  fconn = cell(numsub(groupcnt),1);
  for i = 1:numsub(groupcnt)
    if isempty(parent_folder{1}{groupcnt})
      secondpart = subjects{groupcnt}{i}(1:2);
      if str2double(secondpart) > 96
        pfolder = ['19' secondpart];
      else
        pfolder = ['20' secondpart];
      end
    else
      pfolder = parent_folder{1}{groupcnt};
    end
    disp('--------------------------------------------------------------');
    fprintf('Processing subject: %s \n', subjects{groupcnt}{i});
    subdir = fullfile(server_path, pfolder, subjects{groupcnt}{i}, ...
      'fmri', 'stats_spm8', stats_folder);
    fconn{i} = subject_fconn_task_ROI(roi_folder, roi_list, subdir, ...
      filt, num_trunc, TR, 1, global_ts);
  end
  %-Get groupstats for ROI functional connectivity
  results_fconn_task_ROI.subjectdata{groupcnt} = fconn;
  results_fconn_task_ROI.withingroup{groupcnt} = group_fconn_task_ROI(fconn);
end

numrois = fconn{1}.nrois;
numsess = length(results_fconn_task_ROI.withingroup{1});
numcond = zeros(numsess,1);
for i = 1:numsess
  numcond(i) = length(results_fconn_task_ROI.withingroup{1}{i});
end
  
if numgroup == 2
  disp('----------------------------------------------------------------');
  disp('Do comparisons between two groups');
  numroipair = size(results_fconn_task_ROI.withingroup{1}{1}{1}.data, 2);
  for sesscnt = 1:numsess
    for condcnt = 1:numcond(sesscnt)
      p_ttest = zeros(numroipair,1);
      p_ranksum = zeros(numroipair,1);
      group1 = results_fconn_task_ROI.withingroup{1}{sesscnt}{condcnt}.data;
      group2 = results_fconn_task_ROI.withingroup{2}{sesscnt}{condcnt}.data;
      for ipair = 1:numroipair
        % t test of two independent samples
        [h, p] = ttest2(group1(:,ipair), group2(:,ipair));
        p_ttest(ipair) = p;
        % wilcox rank-sum test of two independent samples
        p = ranksum(group1(:,ipair), group2(:,ipair));
        p_ranksum(ipair) = p;
      end
      results_fconn_task_ROI.betweengroup{sesscnt}{condcnt}.ttest = ...
        tril2full(p_ttest, numrois);
      results_fconn_task_ROI.betweengroup{sesscnt}{condcnt}.ranksum = ...
        tril2full(p_ranksum, numrois);
    end
  end
end

%-Write data to .txt files
disp('------------------------------------------------------------------');
fprintf('Writing data into .txt files ......\n');
rois  = results_fconn_task_ROI.subjectdata{1}{1}.roi_name;
nrois = length(rois);
colname = cell((nrois)*(nrois-1)/2+1,1);
colname{1} = '(Subjects)\(ROI Pairs)';
count = 2;
for i = 1:nrois
  for j = i+1:nrois
    colname{count} = ['(', rois{i}, ')', 'Vs', '(', rois{j}, ')'];
    count = count + 1;
  end
end

conds = results_fconn_task_ROI.subjectdata{1}{1}.event_name{1};
for gcnt = 1:numgroup
  rowname = subjects{gcnt};
  for scnt = 1:numsess
    for ccnt = 1:numcond(scnt)
      event = conds{ccnt};
      fname = ['group', num2str(gcnt), '_session', num2str(scnt), ...
               '_', event, '.txt'];
      dmatrix = results_fconn_task_ROI.withingroup{gcnt}{scnt}{ccnt}.data;
      write_txt(fname, colname, rowname, dmatrix);
    end
  end
end
        
%-Save results
cd(currentdir);
disp('------------------------------------------------------------------');
disp('Saving results');
save results_fconn_task_ROI.mat results_fconn_task_ROI;

fprintf('Changing back to the directory: %s \n', currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('ROI Functional Connectivity finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');

diary off;
clear all;
close all;

end