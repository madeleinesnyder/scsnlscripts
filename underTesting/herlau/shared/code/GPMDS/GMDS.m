function [p,stats] = GMDS(p,opts)
% rng(1);  
if nargin < 1,  
    test(); return; 
end 
for j=1:length(p.s),
    p.s(j).s2x = mk_s2x(size(p.s(j).y,2), size(p.s(j).Phi,2)); 
end
cdir = fileparts(mfilename('fullpath'));
addpath(fullfile(cdir, '../../matlabfrag')); % only for subtightplot. 

if nargin < 2, opts = struct(); end
close all;
dopts.T = 500;
dopts.plot = true; 
dopts.TB = 1; 
dopts.TPhi = 1;  
dopts.Thypers = 1;   
dopts.Ts = 1; 
dopts.Tsf = 0; 
dopts.TC = 1; 
dopts.TU = 1; 
dopts.stat_iter = 1; 
dopts.debug = false;

opts = setstructfields(dopts,opts);
if isfield(p, 'stats')
    stats = p.stats;
end
stats.ssf = {'s','yhat','B','Phi','U'};
stats.sf = {'C','Pw_sigma','Pe_sigf','Pe_sign','Pe_l', 'Pe_mu','C_eig','Pw_std', 'lp'}; 
stats.ssf_rem = {'s', 'yhat','B','Phi'};

p.debug = opts.debug;
for k=1:length(p.s)
    p.s(k).debug = p.debug; 
    
    if ~isfield(p.s(k),'y_impute'),
        p.s(k).y_impute = false(1, size(p.s(k).y,2));
    end
end

p = MDS_fward(p);
MDS_logp(p);

lp = MDS_logp(p);  
stps = [0];  
t0 = tic();

[M,~,J] = size(p.C); 
N = length(p.s); 
[L,pp] = size(p.s(1).Phi);
T = size(p.s(1).y,2);
fprintf('M=%i,T=%i,N=%i; J=%i,p=%i,L=%i\n',M,T,N,J,pp,L);
if p.debug,
    disp('Debug mode active');
end 
for i=1:opts.T,     
    if p.debug, MDS_selfcheck(p); end 
     
    p = MDS_fward(p);    
    if p.debug, MDS_selfcheck(p); end 
    
    p = MDS_fward(p);    
    for g=1:opts.Ts, p = MDS_sample_s(p); if p.debug, MDS_selfcheck(p); end, end
    if opts.plot,  lp(end+1) = MDS_logp(p); stps(end+1) = 2; end    
    
    for g=1:opts.Tsf, p = MDS_sample_sf(p); if p.debug, MDS_selfcheck(p); end, end
    if opts.plot,  lp(end+1) = MDS_logp(p); stps(end+1) = 2; end
        
    for g=1:opts.TC, p = MDS_sample_C(p); if p.debug, MDS_selfcheck(p);  end, end  
    if max(max(abs(p.C(:,:,1)))) > max(abs(p.PC(1).IV(:))),
        p.C
        disp('failed C constraint. cant happen. ');
        disp('Delta is: ')
        deltaC= max(max(abs(p.C(:,:,1)))) - max(abs(p.PC(1).IV(:)));
        disp(deltaC)
        if abs(deltaC) > 0.1,
            disp('Delta C is large. Time to crash');
            assert(false);
        else
            disp('Delta C is small; likely rounding issue. Continuing...');
        end
    end
    
    if max(abs(p.C(:))) > 10,
        disp('Failure in sampling: C has diverged. abort-abort.');
        stats.SAMPLING_FAILED = 'THIS SAMPLER HAS FAILED DUE TO DIVERGING C';
        p.SAMPLING_FAILED = 'THIS SAMPLER HAS FAILED DUE TO DIVERGING C';        
        break;
    end
    %% check that C matches limitations... Perha ps best to add this as a function in NGAM?
    for j=1:length(p.PC),
        if ~p.PC(j).satisfy_constraints(vec(p.C(:,:,j))),
            disp('FAILURE, FAILED TO SATISFY CONTRAINTS...');
            p.SAMPLING_FAILED = 'THIS SAMPLER failed to satisfy contraints.';   
            stats.SAMPLING_FAILED = p.SAMPLING_FAILED;
            p.C
            break
        end
    end
    p.lp = MDS_logp(p); 
    for g=1:opts.TU, p = MDS_sample_U(p); if p.debug, MDS_selfcheck(p); end, end
    
    for g=1:opts.TPhi, p = sample_Phi(p); if p.debug, MDS_selfcheck(p); end, end
    
    for g=1:opts.TB*1, p = MDS_sample_b(p); if p.debug, MDS_selfcheck(p); end, end    
      
    for g=1:opts.Thypers, p = MDS_sample_hypers(p,opts); if p.debug, MDS_selfcheck(p); end, end
    %%
    p = sample_y(p);
    
    p.Pw_sigma = sqrt(diag(p.Pw(1).S));
    
    
    p.Pe_sigf = [p.Pe.sigf]; 
    p.Pe_sign = [p.Pe.sign];
    p.Pe_l = [p.Pe.l];
    p.Pe_mu = [p.Pe.mu];
    p.Pe_mu = p.Pe_mu(1,:); 
    p.Pw_std = std([p.s.w]'); 
    %% record eigenspectrum.
    p.C_eig  = [];     
    for j=1:size(p.C,3), 
        p.C_eig = [p.C_eig ; eig(p.C(:,:,j))];
    end
    %%
    if mod(i, opts.stat_iter) == 0, 
        stats = rstats(p,stats,stats.sf);    
        %%
        if ~isfield(stats,'s'),
            for k=1:length(p.s),
                stats.s(k) = struct();
            end
        end
        sa = num2cell(stats.s);
        for k=1:length(p.s),
            if ~isfield(stats,'s') || length(stats.s) < k, ss = struct(); else ss = stats.s(k); end
            sa{k} = rstats(p.s(k),ss,stats.ssf);
        end     
       stats.s = [sa{:}];        
%%
        if mod(i,5) == 0 && opts.plot,
            pstats(p,stats); pause(0.1); drawnow;
        end
    end
    Cbeta_mean = sqrt(p.PC(1).th_sigma_GAM.beta); 
    wbeta_mean = sqrt(p.Pw(1).th_sigma_GAM.beta); 
    tpi = toc(t0)/i; trem = (opts.T-i)*tpi; 
    fprintf('%i,S=%i, %g/%g rem. C.beta=%4.4f, w.beta=%4.4f\n',i, length(p.s), tpi, trem,Cbeta_mean,wbeta_mean);        
end
if opts.plot, pstats(p,stats); end
%%
if opts.T >= opts.stat_iter && i == opts.T
    stats = rstats(p,stats,stats.sf,stats.ssf_rem);   
    ss = stats.s;
    stats = rmfield(stats, 's');
    for k=1:length(p.s),
        stats.s(k) = rstats(p.s(k),ss(k),stats.ssf,stats.ssf_rem);
        p.s(k).s2x = [];        
    end    
end
end
function stats = rstats(p,stats,sf,sfrem)
DOREM = nargin > 3; 
for j=1:length(sf), 
    s = sf{j};  
    M = p.(s);
    n = numel(M);
    if ~isfield(stats,s),
        i = 1; 
    else
        i = size(stats.(s),2)+1;
    end
    if DOREM,
        x = stats.(s);    
        x = x(:,ceil(end/2):end);
        stats.([s,'_m']) = reshape(mean(x,2),size(M));
        stats.([s,'_std']) = reshape(std(x,1,2),size(M));                
        if any(strcmp(s,sfrem)),
            stats = rmfield(stats,s);
        end        
    else
        stats.(s)(:,i) = M(:);
    end
    
end
stats.sf = sf;
end