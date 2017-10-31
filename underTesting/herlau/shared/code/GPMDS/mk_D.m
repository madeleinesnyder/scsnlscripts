function P = mk_D(M,sigma0,sig_scale,alpha,beta)
if nargin < 3, 
        sig_scale = 1; 

end
if nargin < 4, 
    alpha = 1; beta = 1; 
end
if nargin < 2, sigma0 = 1; end
p.sig_scale = sig_scale;

P.D = eye(M)*(sigma0*sig_scale)^2;
P.S = P.D;
P.S_inv = diag(1 ./ diag(P.D));
P.mu = zeros(M,1);

P.alpha = alpha;
P.beta = beta;
P.type = 4; 
P.name = 'Diagional MVN D=Sigma';
end