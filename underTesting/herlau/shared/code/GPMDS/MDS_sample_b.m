function p = sample_b(p)
if nargin < 1, GMDS(); return; end

for k=1:length(p.s),
    M = size(p.C,1); 
    for m=1:M,
        p.s(k)  = dsample_b(p,k,m);        
    end    
end
 
end
function ps = dsample_b(p,k,m)
%%
ps = p.s(k);

y = ps.y(m,:)-p.Pe(m).mu';
M = size(ps.y,1);

x = permute( ps.x(:,m,:), [1, 3, 2]);
Phi = ps.Phi;

SB_inv = p.PB.S_inv * eye(M);
muB = p.PB.mu + zeros(M,1);

%%
HH = p.Pe(m).L'\x' * Phi';
S_inv = HH' * HH + SB_inv;

%%
q = Phi * x * p.Pe(m).S_inv * y' + SB_inv * muB;

if p.debug,     
    b = ps.B(m,:)';
    de = ps.e(m,:)' - p.Pe(m).mu; 
    s0 = -1/2 * de' * p.Pe(m).S_inv * de -1/2 * (b - muB)' * SB_inv * (b - muB);    

    C = y * p.Pe(m).S_inv * y' + muB' * SB_inv * muB;
    
    s1 = -1/2 * b' * S_inv * b + q' * b -1/2 * C; %- p.NWB.mu)' * p.NWB.S_inv * (bm-p.NWB.mu);
    assert(eqc(s0,s1)); 
    po = p;  
    bo = b';
end
%%
    
ps.B(m,:) = imvnrnd(S_inv\q,S_inv);
ps.yhat = permute( sum(bsxfun(@times, ps.B', mmx('mult',ps.Phi,ps.x)),1),[2,3,1]  ); 
ps.e = ps.y - ps.yhat;
if p.debug,
    %%
    p2 = p; 
    p2.s(k) = ps; 
    MDS_selfcheck(p2);
    MDS_dlp_test(p2,po,ps.B(m,:)',bo',S_inv\q,S_inv);    
end 

end

function ps = dsample_b_old(p,k,m)
ps = p.s(k);
[M,T] = size(ps.y);
ep = 1*((1:M)' == m);
N = eye(M) - ep * ep';

S_inv = p.PB.S_inv;
q = p.PB.S_inv * p.PB.mu; 
K = -0.5 * p.PB.mu' * p.PB.S_inv * p.PB.mu;
Xi = zeros(M,T);
yh = bsxfun(@minus, p.s(k).y, p.Pe.mu);
for t=1:T, % Xi and yhat are the same....
    phi = ps.Phi * ps.x(:,:,t);
    Xi(:,t) = sum(bsxfun(@times, ps.B', ps.Phi * ps.x(:,:,t) ),1)';
    Xi(:,t) = N * Xi(:,t);
    S_inv = S_inv + phi(:,m) * phi(:,m)' * p.Pe.S_inv(m,m);
    
    q = q + -phi(:,m) * p.Pe.S_inv(m,:) * Xi(:,t) + phi(:,m) * p.Pe.S_inv(m,:) * yh(:,t);
    if p.debug,
        K = K - 0.5 * Xi(:,t)' * p.Pe.S_inv * Xi(:,t) - 0.5 * yh(:,t)' * p.Pe.S_inv * yh(:,t) + Xi(:,t)'*p.Pe.S_inv * yh(:,t);
    end
end
if p.debug,     
    bm = ps.B(m,:)';     
    s0 = -1/2 * (bm - p.NWB.mu)' * p.NWB.S_inv * (bm-p.NWB.mu);
    for t=1:T,
        dw = zeros(M,1);
        for j=1:M,
            dw(j,1) = yh(j,t) - ps.B(j,:)* (ps.Phi * ps.x(:,j,t) );
        end
        s0 = s0 - 1/2 * dw' * p.NWe.S_inv * dw;
    end
     s1 = -0.5 * (bm - S_inv\q)' * S_inv * (bm - S_inv\q)  + 0.5 * q' * (S_inv\q) + K; 
     if abs(s1 - s0)/abs(s0) > 10^-5,
         abs(s1 - s0) 
         assert(false);
     end     
     po = p; 
     bo = ps.B(m,:);
end
ps.B(m,:) = imvnrnd(S_inv\q,S_inv);
% ps.B(m,:) = mvnrnd(q,inv(S_inv));

% fix the errors.
%%
ps.yhat = permute( sum(bsxfun(@times, ps.B', mmx('mult',ps.Phi,ps.x)),1),[2,3,1]  ); 
ps.e = ps.y - ps.yhat;


if p.debug,
    %%
    p2 = p; 
    p2.s(k) = ps; 
    MDS_selfcheck(p2);
    MDS_dlp_test(p2,po,ps.B(m,:)',bo',S_inv\q,S_inv);    
end
end