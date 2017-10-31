function [S_inv_chol,q,S_inv,mu,sUU] = MDS_mk_S_matrix(p,k)
if nargin < 1, test(); return; end
C = p.C;
[M,~,J] = size(C);
T = size(p.s(k).v,2);
U = p.s(k).U;
v = p.s(k).v;
%%
dudt = p.s(k).du; 
Sw_inv = p.Pw.S_inv; 
if p.debug,  
    Ct = zeros(M,M,T-1);
    s = p.s(k).s;
    for t=1:T-1,
        Ct(:,:,t) = eye(M); 
        for ii=1:J,
            Ct(:,:,t) = Ct(:,:,t) + C(:,:,ii)*v(ii,t+1)*p.dT;
        end
    end
    Ct = cat(3,zeros(M),Ct);
    s0 = 0;
    dd0 = [];
    for t=1:T,
        if t == 1, dd0(:,t) = s(:,t) - U * dudt(:,t);
        else dd0(:,t) = s(:,t) - Ct(:,:,t) * s(:,t-1) - U * dudt(:,t); end
        s0 = s0 -1/2 * trace( dd0(:,t)' * Sw_inv * dd0(:,t));
    end
end
Rb = lagm(T);
L = chol(Sw_inv);
if ~(isfield(p.s(k),'c') && isfield(p.s(k).c,'Spr')),   
   pp = M; q = M; 
   r = 1; ss = T;
   p.s(k).c.Spr = PS(pp,r);
   p.s(k).c.Sqs = PS(q,ss);  
end
Spr = p.s(k).c.Spr;
Sqs = p.s(k).c.Sqs;

LT = L * Spr';
S_inv_chol = kron(eye(T),LT)  - kron( Rb',LT); 
Ct12 = bsxfun(@plus, C * p.dT, 0);    
for ii=1:J,
    Rbv = Rb; 
    Rbv(Rb ~= 0) = v(ii,2:end);
    S_inv_chol = S_inv_chol- kron(Rbv', LT *( Ct12(:,:,ii)));    
end
[I,~] = find(Sqs);
%%
S_inv_chol = S_inv_chol(:,I);
q = zeros(M*T,1);
m = Rb * dudt(:,:)' * U' * Sw_inv;
for j=1:J, 
    Vj =  diag(v(j,:));
    m = m + Rb * Vj * dudt(:,:)' * U' * Sw_inv  * Ct12(:,:,j);
end
q = q - m(:);
m = dudt(:,:)' * U' * Sw_inv ;
q = q + m(:);     

if nargout > 2 || p.debug, 
    sUU = -1/2 * trace( (U * dudt)' * Sw_inv * (U * dudt));
end

if p.debug,
    %%
    vs = vec(s');
    S_inv_chol2 = kron(L, eye(T))  - kron(L, Rb'); 
    UU = U * dudt;
    q2 = UU' * Sw_inv - Rb * UU' * Sw_inv;        
    for j=1:J,
        V = diag(v(j,:));
        S_inv_chol2 = S_inv_chol2 - p.dT * kron(L * p.C(:,:,j), V' * Rb');
        q2 = q2 - p.dT * Rb * V * UU' * Sw_inv * p.C(:,:,j);        
    end
    q2 = vec(q2); 
    s3 = -1/2* vs' * (S_inv_chol2' * S_inv_chol2) * vs + q2' * vs - 1/2 * trace(UU' * Sw_inv * UU);
     
    s1q = q' *vs;
    SS_inv12 = S_inv_chol' * S_inv_chol;    
    s1 = -1/2 * vs' * SS_inv12 * vs + s1q + sUU;
%     [s0, s1, s3] - s0
    if ~eqc(s0,s1) || ~eqc(s0,s3),
         [s0, s1, s3] - s0
        assert(false);
%        34  
    end
end

if nargout > 2, 
    S_inv = S_inv_chol'*S_inv_chol;
    mu = S_inv\q; 
end
end