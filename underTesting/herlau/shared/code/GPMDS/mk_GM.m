function GM = mk_GM(sz,alpha,beta)
if nargin < 1, test(); return; end
if nargin < 2, 
    alpha = 2; % promote sparsity. higher values, more sparse. 
    beta = 1; 
end

GM = struct();
GM.alpha = alpha;
GM.beta = beta; 

GM.S_inv = gamrnd(alpha, 1./beta, sz); % comes from a gamma distribution. Larger is higher values. 
GM.mu = zeros(size(GM.S_inv));
% GM.S_inv = 1 ./ GM.S;
GM.type = 3; 
GM.name = 'Matrix-Gamma';
end