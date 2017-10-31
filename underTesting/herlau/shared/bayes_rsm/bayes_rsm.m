function bayes_rsm()
rng(42); addpath('../code/GPMDS/dists'); addpath('../code/GPMDS/');
clc; close all;
disp('Generating data...');
p = HPC_soc_cog();
opts.T = 4000; 
p =  rsm_mcmc(p,opts); 
end
function [p, stats] = rsm_mcmc(p,opts)
disp('Initializing sampler..');
p = init_rsm(p); 
YH = mk_XHAT(p);
if false
    v = ppar2vec(p);
    mu1 = reshape(YH * v, [size(p.Y)]);
    mu2 = mu2 + p.beta0 + p.alpha(1) * p.Y * p.L+ p.alpha(2) * p.Y * p.L * p.L;
    for k=1:p.K, mu2 = mu2 + bsxfun(@times, (p.betak(k) + p.p(:,k) ), p.X(:,:,k)); end
    for j=1:p.J, mu2 = mu2 + p.s(j) * p.x(:,:,j); end
    mean(abs(mu1(:)-mu2(:)))
end
%%
tic()
YH_sq = (YH' * YH); 
vY_YH = vec(p.Y)' * YH; 
stats = struct();
disp('Running sampler...'); t0 = tic();
sf = {'betak','beta0','alpha','p','s'}; 
for t=1:opts.T
    %% Sample theta here.
    [~,S_inv_beta] = ppar2vec(p);
    S_inv = YH_sq * p.Py.psigma.S_inv + diag(S_inv_beta);
    q = (vY_YH * p.Py.psigma.S_inv)';
    mu = S_inv\q; 
    beta2 = imvnrnd(mu,S_inv);
    p = vec2ppar(p,beta2);
    %% now sample the different elements (i.e. cauchy distributions).
    p.Palpha1 = p.Palpha1.MCMC(p.alpha(1));
    p.Palpha2 = p.Palpha2.MCMC(p.alpha(2));
    for k=1:p.KS
        p.Ps(k) = p.Ps(k).MCMC(p.s(p.sz == k,:));
    end
    for k=1:p.K
        p.Ppk(k) = p.Ppk(k).MCMC(p.p(:,k));    
    end
    %% record stats.     
    for j=1:length(sf), 
        ss = sf{j}; 
        if ~isfield(stats,ss), stats.(ss) = zeros(numel(p.(ss)),p.T); end
        stats.(ss)(:,t) = p.(ss)(:);
    end    
    fprintf('%i of %i; Trem=%.1f\n',t, opts.T, (opts.T-t)*toc(t0)/t );    
end
toc()
%%
close all;
[H,W] = getHW(length(sf)); 
for k=1:length(sf),
    figure(1); subplot(H,W,k); 
    ss = sf{k}; 
    m = stats.(ss);      
    m = m(1:min(5, size(m,1)), :);
    plot(m(:,10:end)');    
    title(ss);   
    %%
    figure(2); subplot(H,W,k); 
    for j=1:size(m,1),
        [cc,xx] = hist(m(j,:),30);  
        cc = cc / trapz(xx, cc);
        plot(xx',cc','.-'); hold all;
    end
    title(ss);
end
mean(stats.betak(:,end/2:end),2)
end
function YH = mk_XHAT(p)
YH = 0;
for i=1:p.N
    ei = ((1:p.N)' == i)*1;     
    Yi = [p.Y(i,:) * p.L ; p.Y(i,:) * p.L * p.L ; ones(1,p.T)];        
    for k=1:p.K
        Xi = p.X(i,:,k);  
        Yi(3+k,:) = Xi; 
    end    
    for j=1:p.J
        xj = p.x(i,:,j);  
        Yi(3+p.K+j,:) = xj; 
    end
    YH = YH + kron(Yi', ei);    
end
dX = 0; 
for i=1:p.N
    Xk = 0;
    ei = ((1:p.N)' == i)*1;
    Xi = p.X(i,:,k);         
        
    for k=1:p.K
        ek = ((1:p.K)' == k)*1;                 
        Xk = Xk + Xi' * ek'; 
    end   
    dX = dX + kron( Xk, ei * ei');    
end
YH = [YH, dX];
end
function [beta,S_inv_beta] = ppar2vec(p)
beta = [p.alpha ; p.beta0 ; p.betak ; p.s ; vec(p.p) ];
S_inv_s = zeros(p.J,1);
for k=1:p.KS
    S_inv_s = S_inv_s + (p.sz == k) * p.Ps(k).psigma.S_inv;     
end
S_inv_beta = [p.Palpha1.S_inv ; p.Palpha2.S_inv ; p.Pbeta0k.S_inv ; S_inv_s ; vec(repmat(arrayfun(@(pk)pk.psigma.S_inv, p.Ppk), [p.N,1])) ];
end
function p = vec2ppar(p,beta)
p.alpha = beta(1:2)';
p.beta0 = beta(3);
p.betak = beta(4:3+p.K)';
p.s = beta(3+p.K+1:3+p.K+p.J)';
p.p = reshape( beta(3 + p.K+p.J +1:end)', size(p.p)); 
end
function p = init_rsm(p)
%% set up distributions of various types. 
[p.N,p.T,p.K] = size(p.X);
p.J = size(p.x,3);
CAUCHYVAL = 1; 
p.Py = th_NhalfC(CAUCHYVAL);
p.Palpha1 = th_cauchy(1);
p.Palpha2 = th_cauchy(1);
p.Pbeta0k.S_inv = ones(p.K+1,1) * (10^6)^-2;
for k=1:p.K
    p.Ppk(k) = th_NhalfC(CAUCHYVAL);
end 
p.KS = length(unique(p.sz));
for k=1:p.KS
    p.Ps(k) = th_NhalfC(CAUCHYVAL);
end
p.L = eye(p.T); % lag matrix. 
p.L = [p.L(:,2:end), p.L(:,1)*0];
p.L = p.L'; 
end
function p = HPC_soc_cog(N)
T = 200;
J = 10; % number of stimuli. 
if nargin < 1, N = 100; end %Number of subjects. 

xij = zeros(N,T,J);
for i=1:N
    for j=1:J
        dx = zeros(1,T);
        for v=1:2, hh = randi(T-20); dx(hh:hh+20) = 1; end
        xij(i,:,j) = dx;
    end
end
p.x = xij;
z = round(linspace(1,2,J) )'; 
K = length(unique(z));
for k=1:K
    X(:,:,k) = sum(p.x(:,:,z == k),3);     
    p.betak(k,1) = k*2;
end 
p.X = X;
p.beta0 = 0; 
p.sz = z; 
p = init_rsm(p);
p.p = randn(p.N,p.K)/10;
for k=1:p.KS
   ss = randn()/10;
   I = p.sz == k;
   p.s(I,:) = rand(nnz(I),1) * ss;
end
p.alpha = [0 ; 0];
p.Py.sigma = 0.05;  p.Py.S_inv = 1/p.Py.sigma^2;
p = rsm_forward(p);
end
function p = rsm_forward(p) % forward simulation. 
p.Y = zeros(p.N,p.T);
alpha = p.alpha;
p.alpha = alpha*0;
XH = mk_XHAT(p);
beta = ppar2vec(p);
mu = reshape(XH  * beta, size(p.Y));
% do autoregressive stuff. 
p.alpha = alpha;
for t=1:p.T
    yb  = p.Y * p.L;
    ybb = p.Y * p.L * p.L;    
    p.Y(:,t) = yb(:,t)*p.alpha(1) + ybb(:,t) * p.alpha(2) + mu(:,t) + randn(p.N,1)*p.Py.sigma;            
end
end