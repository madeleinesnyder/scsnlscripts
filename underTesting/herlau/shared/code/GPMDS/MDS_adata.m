function [p,yhat] = MDS_adata(p)
if nargin < 1, 
    M = 3; 
    C = zeros(M,M);
    C(2,1) = .5; 
    C(1,2) = -c(2,1);    
    C(3,1) = .9;    
    C = C - eye(M)/4;
    dT = 0.05; 
    T = 500;
    v = ones(1,T);         
    p.C = C;
    p.Pw.D = eye(M)/10^2;
    p.Pw.mu = zeros(M,1); 
    p.Pe.D = eye(M)/20^2;
    p.Pe.mu = zeros(M,1);          
    p.dT = dT;
    p.s(1).v = v;     
    Phi_mu = spm_hrf(p.dT);
    pp = 3;  
    p.s(1).Phi = repmat(Phi_mu',[pp,1]);
    p.s(1).B = eye(M,pp);  
    p.s(1).u = zeros(M,T);
    p.s(1).U = zeros(M);
end
[S_inv_chol,q] = mk_S_matrix(p.C,p.dT,p.Pw.D,p.s(1).v);
[~,T] = size(p.s(1).v);
[M,~,J] = size(p.C);
S_inv = S_inv_chol' * S_inv_chol;
mu = (S_inv) \ q;
r = randn(1,size(S_inv_chol,1))/(S_inv_chol') + mu';
s = reshape(r',[T,M])';
HH = 4; WW = 1;
subplot(HH,WW,1);
plot(s')
%% make Y.
H = p.s(1).B * p.s(1).Phi;
yhat = zeros(M,T);
L = size(H,2);
for t=1:T, 
    Leff = min(L,t);    
    dH = H(:,1:Leff);     
    yhat(:,t) = sum(dH(:,end:-1:1) .* s(:,t-Leff+1:t),2);    
end

subplot(HH,WW,2);
plot(yhat');

h2 = cat(2,H(:,1:min(end,T)),zeros(M,max(0,T-L)));
subplot(HH,WW,3);
plot(h2');
p.s(1).y = yhat + mvnrnd(p.Pe.mu, p.Pe.D,T)';  % D is the covariance matrix. Squared elements. 
subplot(HH,WW,4);
plot(p.s(1).y');

p.s(1).s_true = s; 
p.s(1).s = s;  
end