function [S_inv_chol,q,sUU] = mk_S_matrix(C,dT,D,v,U,u)
p.debug = false;

if nargin < 4,    
    clc
    M = 3; 
    J = 2; 
    T = 240; 
    T = 10;
    
    p.debug = true;
    C = rand(M,M,J);
    dT = 0.1; 
    D = diag( rand(M,1) );    
    v = rand(J,T);
    if J==1,
        v(:) = 1;
    end
end
[M,~,J] = size(C);
T = size(v,2);
if nargin < 6,
    U = diag(rand(M,1))*p.debug;
    u = rand(M,T)*p.debug;
end
if p.debug, 
    s = rand(M,T);    
end


dudt = [u(:,2:end)-u(:,1:end-1),zeros(M,1)];
%%  
 
SD_inv = inv(D * dT);

if p.debug,
    Ct = zeros(M,M,T-1);

    for t=1:T-1,
        Ct(:,:,t) = eye(M);
        for ii=1:J,
            Ct(:,:,t) = Ct(:,:,t) + C(:,:,ii)*v(ii,t)*dT;
        end
    end
    Ct = cat(3,zeros(M),Ct);

    s0 = 0;
    dd0 = [];
    for t=1:T,
        if t == 1, dd0(:,t) = s(:,t) - U * dudt(:,t);
        else dd0(:,t) = s(:,t) - Ct(:,:,t) * s(:,t-1) - U * dudt(:,t); end
        s0 = s0 -1/2 * trace( dd0(:,t)' * SD_inv * dd0(:,t));
    end
end

Rb = circshift(eye(T),-1,1);
Rb(end,1)=0;

%%
L = chol(SD_inv);
pp = M; q = M; 
r = 1; ss = T;
Spr = PS(pp,r);
Sqs = PS(q,ss);
LT = L * Spr';

S_inv_chol = kron(eye(T),LT);
Ct12 = bsxfun(@plus, C * dT, eye(M));    
for ii=1:J
    Rbv = Rb; 
    Rbv(Rb ~= 0) = v(ii,1:end-1);
    S_inv_chol = S_inv_chol - kron( Rbv',LT * Ct12(:,:,ii));        
end     
[I,~] = find(Sqs);
S_inv_chol = S_inv_chol(:,I);    
    
q = zeros(M*T,1);
for j=1:J,
    Vj =  diag(v(j,:));
    m = Vj * Rb * dudt(:,:)' * U' * SD_inv  * Ct12(:,:,j);
    q = q - m(:); 
end
m = dudt(:,:)' * U' * SD_inv ;
q = q + m(:);     

if nargout > 2 || p.debug, 
    sUU = -1/2 * trace( (U * dudt)' * SD_inv * (U * dudt));
end


if p.debug,
    X = s'; vX = X(:);
    s1q = q' *vX;
    SS_inv12 = S_inv_chol' * S_inv_chol;    
    s1 = -1/2 * vX' * SS_inv12 * vX + s1q + sUU;
    
    [s0,s1]-s0
end


end

function S = PS(m,n)
In = eye(n);

for i=1:m,
    ei = (1:m)' == i;
    dm =     kron(ei', kron(In,ei)); 
    if i == 1,
        S = dm;
    else
        S = S + dm;
    end

end
end
