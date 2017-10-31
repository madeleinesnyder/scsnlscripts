function NW = mk_NW(M,kappa,nu,t)
if nargin < 1, test(); return; end
NW = struct();
NW.mu = zeros(M,1);
NW.S_inv = eye(M);
% 2
if nargin < 2, 
    kappa = 1;
end
if nargin < 3, 
    nu = M * 2; % must satisfy nu > 2M-2
    t = 1;
end

NW.kappa = kappa;
NW.nu = nu;

if size(t,1) ~= M, 
    NW.T = eye(M) * t;  
else
    NW.T = t;
end
NW.mu0 = zeros(M,1);

NW.type = 1; 
NW.name = 'Normal-Wishart';
end