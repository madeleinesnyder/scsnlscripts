function [stats,p2,opts] = GPMDS_matlab(p,opts,Bval,Cval)
opts.T = min(opts.T, 2500);
if nargin ==0,
    oddball_GPMDS_convergeA();
    return;
end
if isfield(opts,'firstrun') && ~opts.firstrun,
    %%
    cdir = fileparts(mfilename('fullpath')); 
    addpath(fullfile(cdir, 'mmx'))  
    mmx(1);
    [p2,stats] = GMDS(p,opts);
    stats.Pw_sigma = [];
    return
end
 
USE_C_sigma2 = true;
if nargin < 4, Cval = 30; end
if nargin < 3, Bval = 3; end

SAMPLE_C_hypers = isnan(Cval);
SAMPLE_w_hypers = isnan(Bval);

od.all_NaN = false;
opts = ssfr(od,opts);
    
cms = 1; Cmodes = Cval; 

TT = opts.T;
[M,T] = size(p.s(1).y);
JJ = size(p.C,3); 
for jj=1:JJ,
    opts.PC(jj).type = 7;  %NGAM
    opts.PC(jj).alpha = 1;
    opts.PC(jj).beta = opts.PC(jj).alpha * Cmodes(cms);      
    
    opts.PC(jj).opts = struct();
    if USE_C_sigma2        
        %opts.PC(jj).sigma2_on = true; 
        opts.PC(jj).opts.sigma2_on = true;
        opts.PC(jj).alpha = opts.PC(jj).alpha^2;
        opts.PC(jj).beta  = opts.PC(jj).beta^2;         
        %sqrt(opts.PC(jj).alpha / opts.PC(jj).beta)
    end
    %%
    %Cmean = opts.PC(jj).beta; 
    %Cstd = 8*10^2;
    %alpha = (Cmean / Cstd)^2; 
    %beta = Cmean / Cstd^2; 
    %x = sqrt(gamrnd(alpha,1/beta,[1000,1]));     
    if SAMPLE_C_hypers && ~opts.all_NaN,
        opts.PC(jj).beta = 20^2; 
    end
    
    [alpha,beta] = gammaCenteredAt(opts.PC(jj).beta, 1/8);
    
    opts.PC(jj).opts.Palpha.sample = 0; 
    opts.PC(jj).opts.Pbeta.sample = SAMPLE_C_hypers || opts.all_NaN;
    opts.PC(jj).opts.Pbeta.alpha = alpha;
    opts.PC(jj).opts.Pbeta.beta = beta;     
    opts.PC(jj).opts.name = 'P(C)';
    %%   
    % generate cutoff here.  
    tau_max = 5; % maximum relaxation time.
    C_cut = 1/p.dT * (exp(-p.dT/tau_max) - 1);    
    I = eye(M);
    opts.PC(jj).opts.IV = vec(I) * [-.5, C_cut] + vec(1-I) * [-1, 1] * .5; %.4, .4];
    p.C(:,:,jj) = rand(M)/10 - eye(M)/3;    
end 
%% P(w)
opts.Pw.sigma2_on = true; 
opts.Pw.opts.sigma2_on = true;
opts.Pw.alpha = 1;     

KC = ( 1-(1+p.dT * 2 * C_cut)^2 );
if SAMPLE_w_hypers && ~opts.all_NaN, 
    Bval = 2; 
end
beta_actual2 = (opts.Pw.alpha / (Bval^2 * KC) );
 
opts.Pw.beta = beta_actual2;
opts.Pw.opts.name = 'P(w)';
if SAMPLE_w_hypers || opts.all_NaN,
    %%
    [alpha,beta] = gammaCenteredAt(beta_actual2, 1/16);    
    opts.Pw.opts.Palpha.sample = 0;
    opts.Pw.opts.Pbeta.sample = SAMPLE_w_hypers || opts.all_NaN;
    opts.Pw.opts.Pbeta.alpha = alpha;
    opts.Pw.opts.Pbeta.beta = beta;     
    opts.Pw.beta = opts.Pw.opts.Pbeta.alpha / opts.Pw.opts.Pbeta.beta;        
end 
opts.Pw.sigma_temporal = true;
%% P(e)
opts.Pe.opts.Tsigf = 0;
opts.Pe.sigf = 0;  

opts.PU.U_meanstim = false;  
p = MDS_init(p,opts);
%%
opts.TPhi = 0;
opts.TB = 0;
opts.TU = 1;

%TT = opts.T;
%opts.T = 10;
%opts.plot = true;
%opts.debug = false;

opts.TC = 1;
if ~isfield(opts,'stat_iter')
    opts.stat_iter = 100;
end

%p2 = GMDS(p,opts); 
%opts.TC = 1; 
%opts.T = TT; 
[p2,stats] = GMDS(p,opts);     
% stats.name = ['U=' num2str(u) ', d=' num2str(d), ' st=' num2str(st), ' pwb=' num2str(pwb) ' cms=' num2str(cms), ' Imds=' num2str(Imds)];
%% 
for k=1:length(p2.Pw),
    stats.Pw(k).sigma_m = mean(p2.Pw(k).stats.sigma(round(end/2):end,:),1); 
    stats.Pw(k).sigma_std = std(p2.Pw(k).stats.sigma(round(end/2):end,:),1,1); 
end
for j=1:length(p2.s),
    stats.p.s(j).s = p2.s(j).s;
    stats.p.s(j).v = p2.s(j).v;
    stats.p.s(j).yhat = p2.s(j).yhat;
    stats.p.s(j).y = p2.s(j).y;    
end
for j=1:length(p2.Pe),
    stats.p.Pe(j).mu = p2.Pe(j).mu;
end
    stats.Pw_sigma = [];    
end
function [alpha,beta] = gammaCenteredAt(m,rho)
if nargin < 2
    rho = 1/4; 
end
%%    m = 10; rho = 1/4;         
%    m = 10; rho = 1/4;     
sd = m*rho;
alpha = (m / sd)^2; 
beta = m / sd^2; 
%% gamma with center at m and std of m * rho where e.g. rho = 0.5; 
end

