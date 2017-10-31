function p = ratdat(Ms,J)
%%
if nargin < 1, Ms = 3; end
if nargin < 2, J = 1; end

% close all; clc;
cdir = fullfile( fileparts(mfilename('fullpath')), 'Data');
%% 
ss = dir([cdir,'/*.mat']);
cc = 1;  
p = struct();
plt = nargout > 1; 
for i=1:length(ss),
%     figure(i)
    s = ss(i).name;
%     s
    x = load([cdir, '/' s]);
    y = x.timeseries;
    
%     x
    if isfield(x,'stim_design'),
        v = x.stim_design;
    else
        v = x.Vm';    
    end
    if plt,
    H=5; W = 2; 
    subplot(H,W,cc); cc = cc+1;
    plot(y');     
    title(s)
    subplot(H,W,cc); cc = cc+1;
    imagesc(v);
%     size(y)
    disp(size(y,2));
    end
    
    p.s(i).y = y;
%     p.s(i).v = v*0+; 
%     if J>1,
    p.s(i).v = [v*0+1 ; v]; 
    if J == 1, 
        p.s(i).v = p.s(i).v(1,:);
    end
    T = size(y,2);
    p.s(i).u = zeros(Ms,T,J);
    
    

end
p.roi_names = {'M1','Thalamus','Insula'};
%%

%x = load(fullfile(cdir, '../BDS_Phit'));
%for k=1:length(p.s),
%    p.s(k).Phi = x.Bv';
%end
p.C = -rand(Ms,Ms,J)/10;
p.dT = 0.75;
%%
end