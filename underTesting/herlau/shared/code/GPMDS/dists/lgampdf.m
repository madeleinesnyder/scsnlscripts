% WIKIPEDIA parameterization.
function y = lgampdf(x,alpha,beta)
if nargin < 2, alpha = 1; beta = 3; end
y = log(beta)*alpha - gammaln(alpha) + (alpha-1) * log(x) - beta * x;
end
