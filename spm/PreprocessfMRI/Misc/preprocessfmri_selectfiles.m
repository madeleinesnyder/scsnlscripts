function [InputImgFile, SelectErr] = preprocessfmri_selectfiles(FileDir, PrevPrefix, DataType)

SelectErr = 0;

switch DataType
  case 'img'
    InputImgFile = spm_select('ExtFPList', FileDir, ['^', PrevPrefix, 'I.*\.img']);
  case 'nii'
    InputImgFile = spm_select('ExtFPList', FileDir, ['^', PrevPrefix, 'I.*\.nii']);
    V = spm_vol(InputImgFile);
    nframes = V(1).private.dat.dim(4);
    InputImgFile = spm_select('ExtFPList', FileDir, ['^', PrevPrefix, 'I.*\.nii'], (1:nframes));
    clear V nframes;
end

InputImgFile = deblank(cellstr(InputImgFile));

if isempty(InputImgFile{1})
  SelectErr = 1;
end

end