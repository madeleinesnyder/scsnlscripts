function y = lhalftpdf(sigma,nu,A)
% if nargin < 3, A = 
% implements logarithm of the half student t distribution.
y = log(2) - log(A) - 1/2 * log(nu) - betaln(1/2, nu/2) ...
    + (-(nu+1)/2).*log(1 + 1/nu * (sigma.^2 / A^2) );
end