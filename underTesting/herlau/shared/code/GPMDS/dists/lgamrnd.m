% Numerically stable sampling of gamma distibution for small shape
% parameters where matlab fails. 
function Y = lgamrnd(alpha,beta,varargin)
if nargin < 1,
    dbug = true;
    beta = 2; 
    alpha = 10^-8; 
    alpha = 0.8;  
    n = 100000; 
    Y = exp(lgamrnd(alpha,beta,[1,n]));
    
    x = linspace(alpha/100,alpha/beta * 10);
    
    [y] = histc(Y,x);
    y = y(1:end-1);
    x = x(1:end-1) + (x(2:end) - x(1:end-1))/2;
    y2 = gampdf(x, alpha, 1/beta);
    close all;
    plot(x,y / trapz(x,y)); hold all;
    plot(x, y2 / trapz(x,y2),'k-');
    
    return;
end
[err, sz] = statsizechk(2,alpha,beta,varargin{:});
if err > 0
    error(message('stats:gamrnd:InputSizeMismatch'));
end
assert(all(alpha>0));
alpha = repmat(alpha,sz ./ size(alpha) );
beta = repmat(beta,sz ./ size(beta) );


lambda = 1 ./ alpha  - 1;
w = alpha ./ (exp(1) * (1-alpha)); 
r = 1 ./ ( 1 + w);

if nargin < 3, 
    sz = size(alpha);
end
% I = true(sz);

Z = zeros(sz);
J = alpha > 0.5;
Z(J) = -log(gamrnd(alpha(J),1)) .* alpha(J);

logc =  -gammaln(alpha+1);
i = 0; 
while ~all(Z),
    U = rand(sz);
    J = U < r;   
    
    z = (-log(U ./r )) .* J + ( log(rand(sz)) ./ lambda ) .* (1-J);    
    
    log_haz = logc - z - exp(-z ./ alpha);    
    log_eaz = logc + log(w) + (z >= 0) .* (-z) + (z < 0) .* (log(lambda) + log(w) + lambda .* z);

    J = ~Z & (exp(log_haz - log_eaz ) > rand(sz));
    Z(J) = z(J);
    i = i+1;
    %disp([i, nnz(Z)/numel(J)]);
    if i > 100,
        keyboard;
    end
end
% disp(i);
Y = -log(beta) - Z ./ alpha;% + log(beta); %Y * beta; 
end
% function ggrow(m,)

% end