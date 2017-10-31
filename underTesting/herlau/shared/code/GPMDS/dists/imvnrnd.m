function [r,T] = imvnrnd(mu,inv_sigma,cases)
if nargin < 1, 
    %%    
    close all;
      
    Sigma = [1 0.8 ; .8 2];
    mu = [1 ; 0.5];
    N = 1000000;
%     N = 10;
    rng(1);
    X1 = mvnrnd(mu,Sigma,N);
    rng(1);
    X2 = imvnrnd(mu,inv(Sigma),N);    
        if N < 1000,
        plot(X1(:,1),X1(:,2),'k.'); hold on;
        plot(X2(:,1),X2(:,2),'ro'); hold on;
    end
    cov(X1)-Sigma
    cov(X2)-Sigma
end
if nargin == 4, assert(false); end


if nargin < 2 || isempty(mu) || isempty(inv_sigma)
    error(message('stats:mvnrnd:TooFewInputs'));
elseif ndims(mu) > 2
    error(message('stats:mvnrnd:BadMu'));
elseif ndims(inv_sigma) > 3
    error(message('stats:mvnrnd:BadSigma'));
end

[n,d] = size(mu);
sz = size(inv_sigma);
if sz(1)==1 && sz(2)>1
    % Just the diagonal of Sigma has been passed in.
    sz(1) = sz(2);
    sigmaIsDiag = true;
else
    sigmaIsDiag = false;
end


% Special case: if mu is a column vector, then use sigma to try
% to interpret it as a row vector.
if d == 1 && sz(1) == n
    mu = mu';
    [n,d] = size(mu);
end

% Get size of data.
if nargin < 3 || isempty(cases)
    nocases = true; % cases not supplied
else
    nocases = false; % cases was supplied
    if n == cases
        % mu is ok
    elseif n == 1 % mu is a single row, make cases copies
        n = cases;
        mu = repmat(mu,n,1);
    else
        error(message('stats:mvnrnd:InputSizeMismatchMu'));
    end
end

% Single covariance matrix
if ndims(inv_sigma) == 2
    % Make sure sigma is the right size
    
    if sz(1) ~= sz(2)
        error(message('stats:mvnrnd:BadCovariance2DSize'));
    elseif ~isequal(sz, [d d])
        error(message('stats:mvnrnd:InputSizeMismatch2DSigma'));
    end
    
    if nargin > 3
        assert(false);
        % sigma has already been factored, so use it.
        r = randn(n,size(T,1)) / T + mu;
    elseif sigmaIsDiag
        assert(false);
        % Just the diagonal of sigma has been specified.
        if any(inv_sigma<0)
            error(message('stats:mvnrnd:BadDiagSigma'));
        end
        t = sqrt(inv_sigma);
        if nargout>1
            T = 1./diag(t);
        end
        r = bsxfun(@times,randn(n,d),1./t) + mu;
    else
        % Factor sigma using a function that will perform a Cholesky-like
        % factorization as long as the sigma matrix is positive
        % semi-definite (can have perfect correlation). Cholesky requires a
        % positive definite matrix.  sigma == T'*T
        [T,err] = cholcov(inv_sigma);
        if err ~= 0
            error(message('stats:mvnrnd:BadCovariance2DSymPos'));
        end
        %%
        r = randn(n,size(T,1))/(T') + mu;
%         cov(r)
        
        %%
    end
    
end
end
