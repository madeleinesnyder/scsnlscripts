%%
% It is assumed that L' * L = S (the covariance matrix).
% mu = M x 1
% sigma = 1 x M
% X = N x M
function y = lnormpdf(X,mu,sigma)
[~,M] = size(X);
assert(size(mu,1) == M);
assert(size(sigma,1) == M);
dx = bsxfun(@minus, X, mu'); 
%%
y = -1/2 * log(2*pi) + bsxfun(@plus, -log(sigma)', -bsxfun(@times, dx .^2,  1./ (2 * sigma.^2)' ));
% y2 = log(normpdf(X,mu,sigma));

% y-y2

end