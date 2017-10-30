% This script computes movement statistics
% Statistics: range, standard deviation and root of mean square errors
%
% To run preprocessing:
% start matlab: ml7spm8, in matlab prompt
% >> movement_stats(config_file)
% _________________________________________________________________________
% 2016 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: movement_stats.m 2009-08-06 $
% 07/08/2016: add RMS for translation and rotation -TC
% -------------------------------------------------------------------------

function movement_stats (Config_File)

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
fname = sprintf('movementstats-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp('==================================================================');
fprintf('Movement statistics start at %d/%02d/%02d %02d:%02d:%02d\n',c);
disp('==================================================================');
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

currentdir = pwd;

Config_File = strtrim(Config_File);
if ~ exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off; return;
end

% Make it a function name by dropping '.m'
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

% Read in parameters
server_path      = strtrim(paralist.server_path);
parent_folder    = strtrim(paralist.parent_folder);
subjectlist      = strtrim(paralist.subjectlist);
exp_sesslist     = strtrim(paralist.exp_sesslist);
preprocessfolder = strtrim(paralist.preprocessfolder);
disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

% -------------------------------------------------------------------------
% Read in subjects and sessions
% Get the subjects, sesses in cell array format
parent_folder = ReadList(parent_folder);
subjects      = ReadList(subjectlist);
numsub        = length(subjects);
sesses        = ReadList(exp_sesslist);
numsess       = length(sesses);

% Display the configurations
fmri_path = cell(numsub,1);
if isempty(parent_folder{1})
  for i = 1:numsub
    secondpart = subjects{i}(1:2);
    if str2double(secondpart) > 96
      pfolder = ['19' secondpart];
    else
      pfolder = ['20' secondpart];
    end
    fmri_path{i} = fullfile(server_path,pfolder,subjects{i},'fmri');
  end
else
  for i = 1:numsub
    fmri_path{i} = fullfile(server_path,parent_folder{1},subjects{i}, 'fmri');
  end
end


%--------------------------------------------------------------------------
% Do statistics
% -------------------------------------------------------------------------
cd(currentdir);
mkdir('movementstats');
cd('movementstats');
mstatspath = pwd;
for j = 1:numsess
  outfilename = fullfile(mstatspath, [sesses{j} '.txt']);
  fid = fopen(outfilename, 'w');
  fprintf(fid,['%s\t' ...
               '%s\t%s\t%s\t%s\t%s\t%s\t' ...
               '%s\t%s\t%s\t%s\t%s\t%s\t' ...
               '%s\t%s\t%s\t%s\t%s\t%s\t' ...
               '%s\t%s\t%s\n'], ...
               'Scan_ID', ...
               'Rang_X','Rang_Y','Rang_Z','Rang_Pitch','Range_Roll','Rang_Yaw', ...
               'Std_Dev_X','Std_Dev_Y','Std_Dev_Z','Std_Dev_Pitch','Std_Dev_Roll','Std_Dev_Yaw', ...
               'RMS_X','RMS_Y','RMS_Z','RMS_Pitch','RMS_Roll','RMS_Yaw', ...
               'RMS_Tran', 'RMS_Rot', 'RMS_Total');
             
  sessname = sesses{j};
  for i  = 1:numsub
    session_dir = fullfile(fmri_path{i},sessname);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');                   
    fprintf('Processing: %s \n', session_dir);
    % Check to see if the session direcotry exists
    if ~exist(session_dir, 'dir')
      fprintf('Session directory does not exist: %s\n',session_dir);
      continue;
    end
    cd(session_dir);
    cd(preprocessfolder);
    unix(sprintf('gunzip -fq %s','rp_I*'));
    dname = dir('rp_I*.txt');
    if ~isempty(dname)
      mvmfile = dname.name;
    else
      mvmfile = '';
    end
    if exist(mvmfile, 'file')
      rp_I = load(mvmfile);
    else
      cd(fullfile(session_dir,'unnormalized'));
      unix(sprintf('gunzip -fq %s','rp_I*'));
      dname = dir('rp_I*.txt');
      if ~isempty(dname)
        mvmfile = dname.name;
      else
        mvmfile = '';
      end
      if exist(mvmfile,'file')
        rp_I = load(mvmfile);
      else
        disp('No movement file is found ...');
        diary off; return;
      end
    end
    cd(mstatspath);
    nrow = size(rp_I,1);
    mrange = [];
    st_dev = [];
    rmsq   = [];
    scanID = subjects{i};
    mrange = range(rp_I);
    st_dev = std(rp_I);
    rmsq   = sqrt(sum(rp_I.^2)/nrow);
    %-----------------------
    %- add summary RMS for translation and rotation
    trans_mvmnt = rp_I(:, 1:3);
    trans_mvmnt_one = sqrt(sum(trans_mvmnt.^2, 2));
    
    rot_mvmnt = 65.*rp_I(:, 4:6);
    rot_mvmnt_one = sqrt(sum(rot_mvmnt.^2, 2));
    
    total_mvmnt = [trans_mvmnt, rot_mvmnt];
    total_mvmnt_one = sqrt(sum(total_mvmnt.^2, 2));
    
    rms_trans = sqrt(sum(trans_mvmnt_one.^2)/nrow);
    rms_rot = sqrt(sum(rot_mvmnt_one.^2)/nrow);
    rms_total = sqrt(sum(total_mvmnt_one.^2)/nrow);
    
    rms_trt = [rms_trans, rms_rot, rms_total];
    %-------------------------
    fprintf(fid,['%s\t' ...
                 '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t' ...
                 '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t' ...
                 '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t' ...
                 '%.6f\t%.6f\t%.6f\n'], ...
                  scanID, mrange, st_dev, rmsq, rms_trt);

  end
  fclose(fid);
end

cd(currentdir);
  
c     = fix(clock);
disp('==================================================================');
fprintf('Movement Statistics finishes at %d/%02d/%02d %02d:%02d:%02d\n',c);
disp('==================================================================');
diary off;
clear all;
close all;
end