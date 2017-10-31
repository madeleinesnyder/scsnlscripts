function GM = mk_HT(M,zero_mu,sigma,s0,nu,A)
if nargin < 1, test(); return; end
if nargin < 2,
    zero_mu = true;
end
if nargin < 3,
    sigma = 1; % constant scalar. 
end
if nargin < 4, 
    s0 = 1; 
end
if nargin < 5,    
    A=1;
    nu = 4;  
end
GM.zero_mu = zero_mu;
GM.s0 = s0;
GM.A = A;
GM.nu = nu;
GM.M = M;

GM.mu = zeros(M,1);
GM.sigma = sigma * ones(M,1); % the actual sigma values.

GM.S = GM.s0^2 * GM.sigma.^2;
GM.S_inv = diag(1 ./ GM.S);
GM.S = diag(GM.S);

GM.type = 5;
GM.name = 'Student half-t distribution';
if false, % only relevant for imputation.
    A = GM.nu/2;
    B = GM.A^2 * GM.nu/2;
    GM.alpha = A;
    GM.beta = B;
    GM.tau = gamrnd(A, 1./B, [M,1]);
end
GM.p_mu.mu0 = zeros(M,1);
GM.p_mu.S0 = eye(M)*2;
GM.p_mu.S0_inv = inv(GM.p_mu.S0);
end
% http://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
% (see equation 18).
% http://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture5.pdf
