% Generate global signals
%__________________________________________________________________________
%-SCSNL, Tianwen Chen, 2012-03-27

function preprocessfmri_globalsignal(ConfigFile)

CurrentDir = pwd;

disp('==================================================================');
fprintf('Current directory: %s\n', CurrentDir);
fprintf('Script: %s\n', which('preprocessfmri_globalsignal.m'));
fprintf('Configfile: %s\n', ConfigFile);
disp('***Send error messages to: tianwenc@stanford.edu***');
fprintf('\n');

ConfigFile = strtrim(ConfigFile);
if ~strcmp(ConfigFile(end-1:end), '.m')
  ConfigFile = [ConfigFile, '.m'];
end

if ~exist(fullfile(CurrentDir, ConfigFile), 'file')
  fprintf('Error: cannot find the configuration file ... \n');
  return;
end

ConfigFile = ConfigFile(1:end-2);
eval(ConfigFile);
clear ConfigFile;

SubjectList     = strtrim(paralist.SubjectList);
SessionList     = strtrim(paralist.SessionList);

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
clear paralist;
disp('==================================================================');

%==========================================================================
%-Hard-coded configurations
%-Must notify tianwenc@stanford.edu if you make any changes below
DataType = 'nii';
OutputFolder = 'smoothed_spm8';
ServerPath = '/fs/musk1';
TemplatePath = '/home/fmri/fmrihome/SPM/spm8_scripts/BatchTemplates';
%==========================================================================
if ~exist(TemplatePath, 'dir')
  disp('Error: template folder does not exist!');
  return;
end

Subjects = ReadList(SubjectList);
NumSubj  = length(Subjects);
Sessions = ReadList(SessionList);
NumSess  = length(Sessions);

NumTotalSess = NumSubj*NumSess;
TotalSessionDir = cell(NumTotalSess, 1);

spm('defaults', 'fmri');
spm_jobman('initcfg');
delete(get(0, 'Children'));

SessCnt = 0;
for iSubj = 1:NumSubj
  YearId = ['20', Subjects{iSubj}(1:2)];
  fprintf('Processing subject: %s\n', Subjects{iSubj});
  
  for iSess = 1:NumSess
    SessCnt = SessCnt + 1;
    fprintf('---> session: %s\n', Sessions{iSess});
    
    TotalSessionDir{SessCnt} = fullfile(ServerPath, YearId, Subjects{iSubj}, 'fmri', ...
      Sessions{iSess});
    
    TempDir = fullfile(TotalSessionDir{SessCnt}, 'temp_gs');
    
    UnnormDir = fullfile(TotalSessionDir{SessCnt}, 'unnormalized');
    
    if ~exist(TempDir, 'dir')
      mkdir(TempDir);
    else
      unix(sprintf('/bin/rm -rf %s', fullfile(TempDir, '*')));
    end
    unix(sprintf('cp -af %s %s', fullfile(UnnormDir, ['I*', DataType, '*']), ...
      TempDir));
    unix(sprintf('gunzip -fq %s', fullfile(TempDir, 'I*.gz')));
    
    
    OutputDir = fullfile(ServerPath, YearId, Subjects{iSubj}, 'fmri', ...
      Sessions{iSess}, OutputFolder);
    
    if ~exist(OutputDir, 'dir')
      mkdir(OutputDir);
    end
    
    OutputLog = fullfile(OutputDir, 'log');
    if ~exist(OutputLog, 'dir')
      mkdir(OutputLog);
    end
    
    PrevPrefix = '';
    
    ListFile = dir(fullfile(TempDir, [PrevPrefix, 'I*.gz']));
    if ~isempty(ListFile)
      unix(sprintf('gunzip -fq %s', fullfile(TempDir, [PrevPrefix, 'I*.gz'])));
    else
      [InputImgFile, SelectErr] = preprocessfmri_selectfiles(TempDir, PrevPrefix, DataType);
      if SelectErr == 1
        disp('no scans selected');
        continue;
      end
      preprocessfmri_realign('swar', CurrentDir, TemplatePath, InputImgFile, TempDir)
    end
    
    ListFile = dir(fullfile(OutputDir, ['rp_', PrevPrefix, 'I*.txt*.gz']));
    if ~isempty(ListFile)
      unix(sprintf('gunzip -fq %', fullfile(OutputDir, ['rp_', PrevPrefix, 'I*.txt*.gz'])));
    else
      ListFile = dir(fullfile(OutputDir, ['rp_', PrevPrefix, 'I*.txt']));
      if isempty(ListFile)
        unix(sprintf('cp -af %s %s', fullfile(TempDir, ['rp_', PrevPrefix, 'I*.txt']), OutputDir));
      end
    end
    
    if strcmpi(DataType, 'img')
      P = spm_select('ExtFPList', TempDir, ['^r', PrevPrefix, 'I.*\.img']);
    else
      P = fullfile(TempDir, ['r', PrevPrefix, 'I.nii']);
    end
    VY = spm_vol(P);
    NumScan = length(VY);
    disp('calculating the global signals ...');
    if exist(fullfile(OutputDir, 'VolumRepair_GlobalSignal.txt'), 'file')
      delete(fullfile(OutputDir, 'VolumRepair_GlobalSignal.txt'));
    end
    
    fid = fopen(fullfile(OutputDir, 'VolumRepair_GlobalSignal.txt'), 'w+');
    for iScan = 1:NumScan
      fprintf(fid, '%.4f\n', spm_global(VY(iScan)));
    end
    fclose(fid);
    unix(sprintf('/bin/rm -rf %s', TempDir));
  end
end


delete(get(0, 'Children'));
clear all;
close all;
disp('==================================================================');

end
