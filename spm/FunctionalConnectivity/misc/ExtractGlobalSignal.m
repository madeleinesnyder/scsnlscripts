function [GSignal] = ExtractGlobalSignal(data_type, imagefilter, datadir)

unix(sprintf('gunzip -fq %s', fullfile(datadir, [imagefilter, '*.gz'])));

switch data_type
  case 'img'
    files = spm_select('ExtFPList', datadir, ['^',imagefilter,'.*\.img']);
    nscans       = size(files,1);            
  case 'nii'
    nifti_file = spm_select('ExtFPList', datadir, ['^',imagefilter,'.*\.nii']);
    V       = spm_vol(deblank(nifti_file));
    if length(V) == 4
      nframes = V.private.dat.dim(4);
      files = spm_select('ExtFPList', datadir, ['^',imagefilter,'.*\.nii'],1:nframes);
    else
      files = nifti_file;    
    end
    nscans = size(files, 1);
    clear nifti_file V nframes;
end

VY = spm_vol(files);
GSignal = zeros(nscans, 1);

for iScan = 1:nscans
  GSignal(iScan) = spm_global(VY(iScan));
end

%unix(sprintf('gzip -fq %s', fullfile(datadir, [imagefilter, '*'])));

end