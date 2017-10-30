function preprocessfmri_smooth(WholePipeLine, TemplatePath, InputImgFile, OutputDir, SmoothWidth)

load(fullfile(TemplatePath, 'batch_smooth.mat'));

matlabbatch{1}.spm.spatial.smooth.fwhm = SmoothWidth;
matlabbatch{1}.spm.spatial.smooth.data = InputImgFile;
BatchFile = fullfile(OutputDir, 'log', ['batch_smooth_', WholePipeLine, '.mat']);
save(BatchFile, 'matlabbatch');
spm_jobman('run', BatchFile);
clear matlabbatch;

end