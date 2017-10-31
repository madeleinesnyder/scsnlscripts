function S = K(GP, x, xprime) 
if nargin < 3, xprime = x; end
dd = bsxfun(@minus, x, xprime');
% kernel function.
if ~isprop(GP,'gptype') || GP.gptype == 1, % Gaussian Process
    S = GP.sigf^2 * exp(-dd.^2 / (2*GP.l^2));
else
    % see http://jmlr.csail.mit.edu/proceedings/papers/v1/archambeau07a/archambeau07a.pdf
    % Ornstein-Uhlenbeck process. 
    S =  GP.sigf^2/2 * GP.l * exp(-abs(dd) / (GP.l));
end
if GP.sigf == 0, 
    S(:) =0;
end 
if nargin < 3, S = S + eye(numel(x)) * GP.sign.^2; 
else S = S + (dd==0) * GP.sign^2 ; end 
end