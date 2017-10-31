function [NW,y,smpl,accept] = NW_sample(NW,D,T_hypers)
if nargin < 1,     
    simex(); return;
end
if nargin < 3, T_hypers = 1; end
%%  
% N = size(D,1);
[kappaN,nuN,TN,muN] = upd_pars(NW,D);

I = inv(TN); I = (I+I') /2;
[~,p] = cholcov(I,1);
if p~=0
   I
   assert(false);
end 

S_inv = wishrnd(I,nuN);
if kappaN == inf,
    mu = muN;
else
    mu = imvnrnd(muN,kappaN * S_inv)';
end
NW.mu = mu;
NW.S_inv = S_inv;

    function y = logf(x)
        if NW.kappa == inf,
            kappa = inf;
            nu = exp(x(1));
            t = exp(x(2));
        else
            kappa = exp(x(1));
            nu = exp(x(2));
            t = exp(x(3));            
        end
        if nu <= 2 * M -2, 
            y = -inf; return;
        end
            
            y = sum(x) + marglike(NW.mu0,kappa,nu,t*eye(M),D);
       if isnan(y) || ~isfinite(y),
          assert(false);
       end
    end 
y = [];
smpl = [];
accept = [];
if T_hypers > 0,
%     clc;
    M = size(NW.S_inv,1);
    y = zeros(1,T_hypers);
%     y(1) = marglike(NW.mu0,NW.kappa,NW.nu,NW.T,D);
    
    pars = log([NW.kappa,NW.nu,NW.T(1)]); %zeros(3,T_hypers+1);    
    if NW.kappa == inf,
        pars = pars(2:end);
    end
    proprnd = @(x)x+randn(size(x))*.1;
    doplot = false;
    %%
    [smpl, accept] = mhsample(pars,T_hypers,'logpdf',@logf,'proprnd',proprnd,'symmetric',1);
     
    for i=1:size(smpl,1),
        y(i) = logf(smpl(i,:));
    end
    if doplot,
        close all;
        subplot(2,2,1); plot(smpl)
        title(accept)

        subplot(2,2,2);
        plot(y)
    end
    
    if NW.kappa == inf,
        xx = exp([inf,smpl(end,:)]);
    else
        xx = exp(smpl(end,:));
    end
    
    NW.kappa = xx(1);
    NW.nu = xx(2);
    NW.T = eye(M)*xx(3);        
end
end
function y = lz(kappa,nu,T),
M = size(T,1);
y = lngengamma(M,nu/2)  - logdet(T) * (nu/2);
if kappa ~= inf, y = y - (M/2)*log(kappa); end
end
function y = marglike(mu0,kappa,nu,T,D)
[N,M] = size(D);
NW.kappa = kappa;
NW.nu = nu;
NW.T = T; 
NW.mu0 = mu0;
[kappaN,nuN,TN] = upd_pars(NW,D);
y = lz(kappaN,nuN,TN) - lz(NW.kappa,NW.nu,NW.T)  - log(2*pi) * (N*M/2);
end
function [kappaN,nuN,TN,muN] = upd_pars(NW,D)
N = size(D,1);
nuN = NW.nu + N;
kappaN = NW.kappa + N;

xbar = mean(D,1)';  

dxb = bsxfun(@minus, D, xbar')';
%%
S = dxb * dxb';
% S2 = S*0;
% for i=1:size(D,1),
%     S2 = S2 + (D(i,:)' - xbar) * (D(i,:)' - xbar)';
% end
% S-S2
%%
TN = NW.T + S + (N / (1 + N/NW.kappa) ) * (NW.mu0 - xbar) * (NW.mu0 - xbar)';
muN = (NW.mu0 + N/NW.kappa  * xbar) / (1 + N/NW.kappa);
end


function simex()

close all;
    Sigma0 = [1 0.9 ; 
        0.9 2];
    mu0 = [2 ; 1]*0;
    D = mvnrnd(mu0,Sigma0,100);
    
    plot(D(:,1),D(:,2),'k.');
    axis equal;
    
%     T_hypers = 0;
    NW = mk_NW(2);
%     NW.kappa = 1;
    
    hold on;
    conf = 0.9;
    h = error_ellipse(Sigma0,mu0,conf);
    set(h,'LineWidth',2,'Color','b');
    
    S = 100;
    T_hypers = 10; 
    prs = zeros(0,2+(NW.kappa < inf));
    y = zeros(1,0);
    accept = zeros(1,0);
    for s=1:S,
        [NW,dy,pars,dacc] = NW_sample(NW,D,T_hypers);

        h = error_ellipse(inv(NW.S_inv),NW.mu,conf);
        set(h,'LineWidth',1,'Color','r');    
        pause(0.1); drawnow; 
        prs = [prs ; pars];  
        y = [y ; dy']; 
        accept = [accept ; dacc];
    end
    %%
    figure(2);
    subplot(2,2,1); 
    plot(prs);
    subplot(2,2,2);
    plot(y);
    
    subplot(2,2,3);
    plot(accept);
    
end