%%
% It is assumed that L' * L = S (the covariance matrix).
function y = lmvnpdf_chol(X,mu,L)
[N,M] = size(X);
dx = bsxfun(@minus, X', mu);
assert(size(mu,1) == M);
% L' * L = S
% S_inv = inv(L) * inv(L')
% L' \ dx
dx = L' \ dx;
LDET = sum(log(diag(L)));
y = -N*M/2 * log(2*pi) - LDET*N - 1/2 * trace(dx' * dx);
end