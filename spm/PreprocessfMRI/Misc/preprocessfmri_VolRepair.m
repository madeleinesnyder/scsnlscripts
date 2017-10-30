function [VolFlag] = preprocessfmri_VolRepair(OutputDir, DataType, ImgPrefix)

addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/ArtRepair'));
VolFlag = 0;

nifti4Dto3D(OutputDir, ImgPrefix);

switch DataType
  case 'img'
    imgfiles = spm_select('FPList', OutputDir, ['^', ImgPrefix, 'I.*\.img']);
    realignfile = spm_select('FPList', OutputDir, '^rp_I_0.*\.txt');
  case 'nii'
    imgfiles = spm_select('FPList', OutputDir, ['^', ImgPrefix, 'I.*\.nii']);
    realignfile = spm_select('FPList', OutputDir, 'rp_I.*\.txt');
end

subflag = scsnl_art_global(imgfiles, realignfile, 1, 2, 0);

if subflag == 1
  VolFlag = 1;
else
  nifti3Dto4D (OutputDir, ['v', ImgPrefix])
end

end