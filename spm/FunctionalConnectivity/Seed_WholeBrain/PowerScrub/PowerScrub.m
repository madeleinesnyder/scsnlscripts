function ScrubIndex = PowerScrub(MvmntMtx, ImgFile)

dmvmnt = [zeros(1, size(MvmntMtx, 2)); abs(diff(MvmntMtx))];
dmvmnt(:,4:6) = dmvmnt(:,4:6)*50;
FD = sum(dmvmnt, 2);
BadMvmntVol = find(FD > 0.5);

MaskFile = '/home/fmri/fmrihome/SPM/spm8_scripts/Masks/vbm_grey_mask.nii';
Mask = spm_read_vols(spm_vol(MaskFile));
MaskIndex = find(Mask > 0);

MaskData = spm_read_vols(spm_vol(ImgFile));
MtxDim = size(MaskData);
MaskData = reshape(MaskData, prod(MtxDim(1:3)), MtxDim(4));
MaskData = MaskData(MaskIndex, :); %#ok<*FNDSB>
MaskData = MaskData';
MaskDataVec = MaskData(:);
MaskDataVec(isnan(MaskDataVec)) = [];
MaskDataVec(MaskDataVec == 0) = [];
MaskMode = max(abs(MaskDataVec));

AvgMaskData = mean(MaskData, 1);
OutBrainIndex = find(AvgMaskData < 0.7*mean(AvgMaskData));

MaskData(:, OutBrainIndex) = [];

DVARS = sqrt(mean(diff(MaskData).^2, 2));
DVARS = [0; DVARS];
DVARS = DVARS./MaskMode*1000;
BadDVARSVol = find(DVARS > 5);

ScrubIndex = intersect(BadMvmntVol, BadDVARSVol);
end