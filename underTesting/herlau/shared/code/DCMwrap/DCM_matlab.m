function [stats,DCM] = DCM_matlab(p,opts)
if nargin < 1, 
    % here we go...
    clc;
    addpath('../');
    p = fMRI_data('attention');    
    z = load('../attention/GLM/DCM_working_copy.mat');
    %% fix the second region link; it is bad, bad bad.
    M = zeros(3); M(2,3) = 1; 
    p.opts.PC(3).IV = vec(M) * [-.5, .5];
    
    %%
end
cdir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(cdir,'../../spm12')));

%%
% z = load('../attention/GLM/DCM_working_copy.mat');
% z.DCM;
DCM = struct();
% DCM.options = struct();
if nargin < 2, opts = struct(); end

DCM.options.nonlinear = 0; 
DCM.options.two_state = 0; 
DCM.options.stochastic = 0; 
DCM.options.centre = 0; 
DCM.options.induced = 0; 

DCM.options.random_init = 0;
DCM.options = ssfr(DCM.options,opts);

% DCM.options.nonlinear = 0;
% DCM.options.two_state = 0;
% DCM.options.stochastic = 1;
% DCM.options.centre = 0;
% DCM.options.induced = 0;

%%
ps = p.s(1);
[M,T]  = size(p.s(1).y);
DCM.v = T;
DCM.n = M;
DCM.TE = 0.04; 
DCM.delays = zeros(M,1)+p.dT / 2;
%%
DCM.d = zeros(M,M,0);
%%
% identify actual stimulation patterns; remove constant irrelevant patterns.
ps_u = ps.u;
ps_v = ps.v; 
% remove constant patterns, e.i. either all 1 or 0.
ps_u = ps_u(:,:,~all(sum(ps_u,1)>0,2));
ps_u = ps_u(:,:,any(sum(ps_u,1),2));

if size(ps_u,3) == 0,
%     ps_u = ones(M,T,1);
end 
J = ~all(ps_v,2);  
if ~isfield(p, 'opts') || ~isfield(p.opts, 'PC'),
    disp('resetting IV for DCM...');
    M = size(p.s(1).y,1);
    for j=1:size(p.s(1).v,1),          
        p.opts.PC(j).opts.IV = ones(M*M,1) * [-1,1];
    end
end
% else
PC = p.opts.PC(J);    
% end


ps_v = ps.v(J,:); %sum(ps_v,
%%
% constant stims have now been removed. 
% identify unique stimulation patterns. 
SPAT = [squeeze(sum(ps_u,1)) ; ps_v];
USPAT = unique(SPAT,'rows');
if size(USPAT,1) == 3, 
    USPAT = USPAT([3,2,1],:);
end
DCM.U.name = cell(1,M);    
for j=1:size(USPAT,1),
    % detect which it matches. 
    % also make entry in U corresponding to this pattern.
%     USPAT(j,:)
    IV = all(bsxfun(@minus, ps_v, USPAT(j,:))==0,2);
    if size(ps_u,3) > 0, 
        IU = all(bsxfun(@minus, squeeze(sum(ps_u,1)>0), USPAT(j,:))==0,2);
    else
        IU = [];
    end
    if ~any(IV),
        DCM.b(:,:,j)= zeros(M);
    else
        % figure out sparsity pattern here. 
        DCM.b(:,:,j) = IV2sparse(PC(IV).opts.IV,M);        
    end
    if ~any(IU), 
        DCM.c(:,j) = zeros(M,1);
    else
        DCM.c(:,j) = sum(ps_u(:,:,IU),2)>0;
    end
    %% make U matrix thingy.
    SCALE = 16; % mystery scale!
    DCM.U.dt = p.dT/SCALE;    
    DCM.U.name{j} = sprintf('Stim %i', j);
    %%
    DCM.U.u(:,j) = zeros(T*SCALE,1);
    for t=1:T*SCALE,
            dx = min(T,-1+ceil( (t)/SCALE));
            dx = max(1,dx);
            DCM.U.u(t,j) = USPAT(j,dx);
    end
end
if isempty(J), 
    DCM.a = zeros(M); 
else
    DCM.a = IV2sparse(p.opts.PC(~J).opts.IV,M);
end

%%
if ~isfield(DCM,'b'),
    DCM.b = zeros(M,M,0);
end
if ~isfield(DCM,'c'),
    DCM.c = zeros(M,0);
end
% maky Y matrix.
DCM.Y.dt = p.dT;
DCM.Y.y = ps.y';
% DCM.Y.y = z.DCM.Y.y;
DCM.Y.name = cell(1,M);
for j=1:M,
    DCM.Y.name{j} = sprintf('region %i',j);
end
%% 
%% matlab crapper.
disp('begin dcm estimation');
% DCM4 = DCM;
disp('a, b, c, d');
DCM.a
DCM.b
DCM.c 
DCM.d

[DCM] = spm_dcm_estimate_tue(DCM); 
stats = struct();
stats.C_m = cat(3,DCM.Ep.A, DCM.Ep.B);
if ~any(DCM.b(:)), 
    stats.C_m = stats.C_m(:,:,1);
end    
stats.F = DCM.F;
stats.qU.x = DCM.qU.x{1};
stats.s.s_m = DCM.qU.x{1}(1:M,:);

return;
%     345;
% else
%     234
% end

% DCM4.b = zeros(M,M,0);
% DCM4.c = zeros(M,0); 
% DCM4.opts.random_init = true;
return;
%%
DCM.b = 1; 
DCM.a = 1; 
DCM.c = 1; 
%%
[dc2] = spm_dcm_estimate_tue(DCM); 
%%
dc1 = spm_dcm_estimate_tue(z.DCM); 
dc3 = spm_dcm_estimate_tue(DCM4); 
%%
close all;
subplot(2,2,1);
imagesc([dc1.Ep.A, dc2.Ep.A,  dc3.Ep.A])
subplot(2,2,2);
cc{1} = flatmat(dc1.Ep.B,.1);
cc{2} = flatmat(dc2.Ep.B,.1);
cc{3} = flatmat(dc3.Ep.B,.1);
imagesc(flatmat(cc));%, dc2.Ep.B,  dc3.Ep.B]))
subplot(2,2,3);
cc{1} = flatmat(dc1.Ep.C,.1);
cc{2} = flatmat(dc2.Ep.C,.1);
cc{3} = flatmat(dc3.Ep.C,.1);
imagesc(flatmat(cc));%, dc2.Ep.B,  dc3.Ep.B]))


%%
dc1 = spm_dcm_estimate(DCM);
dc2 = spm_dcm_estimate(z.DCM);

%%
DCM2 = z.DCM;
DCM2.Y.y = zscore(DCM2.Y.y);
dc3 = spm_dcm_estimate(DCM2);
%%



%%



% imagesc([dc1.Ep., dc2.Ep.A,  dc3.Ep.A])


%%
UM = squeeze(sum(ps.u,1))' >0;
UR = unique([UM ;  ps.v], 'rows');
UR = UR(end:-1:1,:);
if ~all(UR(1,:) == 1)
    assert(false); % no B matrix.
end
DCM.U = struct();
    
for j=1:size(UR,1),
    u = UR(j,:);    
    for i=1:size(ps.u,3),
        
        if all(sum(ps.u(:,:,i),1)>0 == u),
            disp('stim found!');
            % make U struct, somewhere. 
            
        end
    end
end
    

% now fix the matrices. 


% DCM
%%
% DCM = z.DCM;
% p.s(1).IV
DCM.Y = rmfield(DCM.Y,'X0');
%%




end
function b = IV2sparse(IV,M),
b = ones(M);
if size(IV,1) >= 0, 
    b = reshape(IV(:,1) ~= IV(:,2),[M,M]);
    
%     b = zeros(M);
end
% M = sqrt(size(IV,1));

end