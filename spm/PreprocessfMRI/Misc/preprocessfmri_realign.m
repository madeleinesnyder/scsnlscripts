% Preprocess fMRI data - realignment
%__________________________________________________________________________
%-SCSNL, Tianwen Chen, 2011-11-06

function preprocessfmri_realign(WholePipeLine, CurrentDir, TemplatePath, ImgFiles, OutputDir)

% Load default batch
load(fullfile(TemplatePath, 'batch_realign.mat'));

% Update batch parameters
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = {};
matlabbatch{1}.spm.spatial.realign.estwrite.data{1} = ImgFiles;

LogDir = fullfile(OutputDir, 'log');
if ~exist(LogDir, 'dir')
  mkdir(LogDir);
end

% Update and save batch
BatchFile = fullfile(LogDir, 'batch_realign.mat');
save(BatchFile, 'matlabbatch');

% Run batch of realignment
cd(OutputDir);
spm_jobman('run', BatchFile);
clear matlabbatch;
cd(CurrentDir);

PSFile = dir(fullfile(OutputDir, 'spm_*.ps'));
PDFFile = fullfile(LogDir, ['realign_spm8_', WholePipeLine, '.pdf']);
if isempty(PSFile)
  ErrMsg{1} = sprintf('Warning: no reaglignment ps file found');
  disp(ErrMsg{1});
  ErrFlag = 1;
else
  if length(PSFile) > 1
    ErrMsg{1} = sprintf('Warning: multiple ps files found');
    disp(ErrMsg{1});
    ErrFlag = 1;
  else
    unix(sprintf('ps2pdf13 %s %s', ...
      fullfile(OutputDir, PSFile(1).name), PDFFile));
    delete(fullfile(OutputDir, PSFile(1).name));
  end
end

end
