function y = limvnpdf(X,mu,S_inv)
[N,M] = size(X);
dx = bsxfun(@minus, X', mu);
assert(size(mu,1) == M);

y = -N*M/2 * log(2*pi) + 1/2 * logdet(S_inv)*N- 1/2 * trace(dx' * S_inv * dx);
end