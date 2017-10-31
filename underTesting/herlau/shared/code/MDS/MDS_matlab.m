%% MAIN MDS caller. 
% use this script to perform MDS analysis (main analysis).
function [stats,BDS] = MDS_matlab(p,opts)
cdir = fileparts(mfilename('fullpath'));
addpath(fullfile(cdir,'/MDS_dependencies'));
addpath(fullfile(cdir,'../GPMDS'));
% run MDS matlab edition.
if nargin == 0,
    p = fMRI_data('rats');
    opts = struct();
end
if nargin <2, opts = struct(); end
%% 
%% assume same stim design; i.e. flatten stim design here.
for i=1:length(p.s),
    %% 
    ps = p.s(i);
    u = ps.u;
    v = ps.v;    
    u2 = zeros(0,size(p.s(i).u,2));
    for xi=1:size(u,1),
        for xj=1:size(u,3),
            u2(end+1,:) = p.s(i).u(xi,:,xj);
        end
    end
    u = u2(sum(u2,2)~=0,:);
%     u = u2;
        
        
    vv = [v ; u];
    vv = unique(vv,'rows');
    T = size(vv,2);
    vv(sum(vv,2)==0 | sum(vv,2) == T,:) = [];
    vv = [ones(1,T) ; vv];
    p.s(i).v = vv; % flatten the stim design; simplify it (hacky but may work).    
    fprintf('MDS stim design size is J=%i\n', size(vv,1));    
end

od.tol = 10^-6;
od.maxIter = 100;
od.random_init = false;

opts = ssfr(od,opts);
  
[stats, BDS] = main_MDS_estimation_TUHE(p,opts);





end