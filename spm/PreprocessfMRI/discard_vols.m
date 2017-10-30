%-Discard specified volumes from raw data (in unnormalized folder)
%__________________________________________________________________________
%-Stanford SCSNL
%-Tianwen Chen, 2011-10-11

function discard_vols (Config_File)

CurrentDir = pwd;

disp('------------------------------------------------------------------');
fprintf('Current directory: %s\n', CurrentDir);
disp('Running discard_vols.m');
disp('Send error messages to tianwenc@stanford.edu');
disp('------------------------------------------------------------------');

Config_File = strtrim(Config_File);

if ~exist(Config_File,'file')
  fprintf('Error: Cannot find the configuration file ... \n');
  return;
end

% Make it a function name by dropping '.m'
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

ServerPath    = strtrim(paralist.ServerPath);
RawDataFolder = strtrim(paralist.RawDataFolder);
SubjectList   = strtrim(paralist.SubjList);
SessionList   = strtrim(paralist.SessionList);
DiscardVols   = paralist.DiscardVols;

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

Subjects = ReadList(SubjectList);
NumSubj  = length(Subjects);
Sessions = ReadList(SessionList);
NumSess  = length(Sessions);

for iSubj = 1:NumSubj
  YearId = ['20', Subjects{iSubj}(1:2)];
  fprintf('-> Processing subject: %s\n', Subjects{iSubj});
  for iSess = 1:NumSess
    fprintf('---> Processing session: %s\n', Sessions{iSess});
    DataDir = fullfile(ServerPath, YearId, Subjects{iSubj}, 'fmri', Sessions{iSess}, ...
      RawDataFolder);
    cd(DataDir);
    ListFile = dir(fullfile(DataDir, '*.nii*'));
    if isempty(ListFile)
      fprintf('Cannot find the raw data in: %s\n', DataDir);
      continue;
    else
      if length(ListFile) > 1
        disp('Find multiple raw data in: %s\n', DataDir);
        continue;
      else
        nifti4Dto3D(DataDir, []);
        ListFile = dir('*.nii');
        NumFull3D = length(ListFile);
        if ~exist('unused', 'dir')
          unix('mkdir unused');
        end
        for iBegin = 1:DiscardVols(1)
          BeginVolNum = iBegin - 1;
          if BeginVolNum < 10
            FileName = ['I_000', num2str(BeginVolNum), '.nii'];
          else
            if iBegin < 100
              FileName = ['I_00', num2str(BeginVolNum), '.nii'];
            else
              FileName = ['I_0', num2str(BeginVolNum), '.nii'];
            end
          end
          unix(sprintf('mv -f %s unused', FileName));
        end
        
        if length(DiscardVols) == 2
          for iEnd = 1:DiscardVols(2)
            EndVolNum = NumFull3D - iEnd;
            if EndVolNum < 10
              FileName = ['I_000', num2str(EndVolNum), '.nii'];
            else
              if EndVolNum < 100
                FileName = ['I_00', num2str(EndVolNum), '.nii'];
              else
                FileName = ['I_0', num2str(EndVolNum), '.nii'];
              end
            end
            unix(sprintf('mv -f %s unused', FileName));
          end
        end
        nifti3Dto4D(DataDir, []);
      end
    end
  end
  disp('----------------------------------------------------------------'); 
end

disp('Done!');
cd(CurrentDir);

end