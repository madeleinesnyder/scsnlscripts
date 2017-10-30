function xY = ppi_regions(ppi_dir, xY)

%-Extract voxel time series in a volume of interest and find the first
% engenvariate as the physiological part
%-Modified from spm_regions.m
%-Tianwen Chen, 04/15/10

% xY     - VOI structure
%       xY.xyz          - centre of VOI {mm}
%       xY.name         - name of VOI
%       xY.Ic           - contrast used to adjust data (0 - no adjustment)
%       xY.Sess         - session index
%       xY.def          - VOI definition
%       xY.spec         - VOI definition parameters
%       xY.str          - VOI description as a string
%       xY.XYZmm        - Co-ordinates of VOI voxels {mm}
%       xY.y            - [whitened and filtered] voxel-wise data
%       xY.u            - first eigenvariate {scaled - c.f. mean response}
%       xY.v            - first eigenimage
%       xY.s            - eigenvalues
%       xY.X0           - [whitened] confounds (including drift terms)
%
% Y and xY are also saved in VOI_*.mat in the SPM working directory

load SPM.mat;
xY.M = SPM.xVol.M;
XYZ  = SPM.xVol.XYZ;
XYZmm = SPM.xVol.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];
[xY, xY.XYZmm, Q] = spm_ROI(xY, XYZmm);
if isempty(xY.XYZmm)
  error('Empty region.')
end

%-Get raw data, whiten and filter 
y        = spm_get_data(SPM.xY.VY, XYZ(:,Q));
y        = spm_filter(SPM.xX.K, SPM.xX.W*y);

%-Remove null space of contrast
if xY.Ic

    %-Parameter estimates: beta = xX.pKX*xX.K*y
    beta  = spm_get_data(SPM.Vbeta, XYZ(:,Q));

    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    y     = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);

end

%-Confounds
xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);

%-Extract session-specific rows from data and confounds
try
    i     = SPM.Sess(xY.Sess).row;
    y     = y(i,:);
    xY.X0 = xY.X0(i,:);
end

% and add session-specific filter confounds
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
end
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH]; % Compatibility check
end

%-Remove null space of X0
xY.X0     = xY.X0(:,any(xY.X0));


%-Compute regional response in terms of first eigenvariate
[m n]   = size(y);
if m > n
    [v s v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u s u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);

%-Set in structure
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;

%-Save
str = ['VOI_' xY.name];
if isfield(xY,'Sess') && isfield(SPM,'Sess')
    str = sprintf('VOI_%s_%i',xY.name,xY.Sess);
end
if spm_matlab_version_chk('7') >= 0
    save(fullfile(ppi_dir, str),'-V6','Y','xY')
else
    save(fullfile(ppi_dir, str),'Y','xY')
end

fprintf('VOI saved as %s\n', ...
    spm_str_manip(fullfile(ppi_dir, [str '.mat']), 'k55'));

end