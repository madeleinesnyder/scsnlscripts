function p = sample_C(p)
if nargin < 1, MDS(); return; end
[M,~,J] = size(p.C);
for j=1:J, 
    for q=1:M,
            p = dsample_C(p,j,q); 
    end 
end  
end 
function p = dsample_C_GP_defunc(p,j,q)      
Cnj = p.C(:,:,1:end ~= j);
M = size(p.C,1); 
XHTCACHE = struct(); 
A = p.C(:,:,j);
Np = eye(size(A)); Np(q,q) = 0; 
b = p.PC.S_inv * p.PC.mu;    
S_inv = p.PC.S_inv;
if p.debug,
    s0 = -1/2 * (A(q,:)'-p.PC.mu)' * p.PC.S_inv * (A(q,:)'-p.PC.mu);        
    ep = 1*((1:M)' == q);
    P = ep * ep';
    N = eye(M) - P;
    K = p.PC.mu' * p.PC.S_inv * p.PC.mu;
    po = p;
    ap = A(q,:)';    
end

for k=1:length(p.s),
    ps = p.s(k); 
    T = size(ps.y,2);
    
    [MU,S] = GPmarg(p.Pw,ps.w,2:size(ps.w,2));
%     SW_inv = inv(S);
    
    %%
    dm = mtimesx(Cnj,ps.s(:,1:end-1) );
    xh = ps.s(:,2:end) - sum(bsxfun(@times, dm,  permute(ps.v(1:end ~= j, 2:end), [3, 2, 1])),3) ;
    xh = bsxfun(@minus, xh, MU ) -  p.U * p.s(k).v(:,2:end);
    xt = bsxfun(@times, ps.v(j,2:end), ps.s(:,1:end-1));    
    XHTCACHE(k).xh = xh;
    XHTCACHE(k).xt = xt;
    XHTCACHE(k).MU = MU;
    
    %%    
    
    L = cholcov(S);
    dd = xt / L;
    m2 =  dd * dd';   
    bb2 = dd * (L' \ xh(q,:)');
    if p.debug,
       %%
       SW_inv = inv(S);
       m1 = xt * SW_inv * xt';
       bb1 =  xt * SW_inv * xh(q,:)';
           
        assert( mean(abs( m2(:)-m1(:))) < 10^-6); 
        
        assert( mean(abs( bb2(:)-bb1(:))) < 10^-6); 
        
    end
    S_inv = S_inv + m2;
    b = b + bb2;
    %%
    if p.debug && false
        for t=1:size(xh,2),
            s0 = s0 - 1/2 * (xh(:,t) - A * xt(:,t))' * p.Pw.S_inv * (xh(:,t) - A * xt(:,t));
        end
        K = K + trace(xt(:,:,1)' * (N*A)' * p.Pw.S_inv(:,:,1) * (N*A) * xt(:,:,1));            
        K = K + trace( xh(:,:,1)' * p.Pw.S_inv(:,:,1) * xh(:,:,1)) ...
            -2* trace(N * A * xt(:,:,1) *xh(:,:,1)' * p.Pw.S_inv(:,:,1) );
    end
end    
if p.debug && false,

    res2 = -1/2 * (ap - S_inv\b)' * S_inv * (ap - S_inv\b) + 1/2 * (b' * (S_inv\b)) - 1/2 * K;
    if abs(res2-s0) > 10^-10,
        assert(false)
    end 
end 
A(q,:) = imvnrnd(S_inv\b,S_inv);     
p.C(:,:,j) = A;  

for k=1:length(p.s),
    p.s(k).w(:,2:end) = bsxfun(@plus, XHTCACHE(k).MU , XHTCACHE(k).xh - A * XHTCACHE(k).xt);        
end    
if p.debug,
    MDS_selfcheck(p);
    MDS_dlp_test(p,po,A(q,:)',ap,S_inv\b,S_inv)    
end
p.C(:,:,j) = A;   
p = MDS_fward(p); 
end


function p = dsample_C(p,j,q)        
Cnj = p.C(:,:,1:end ~= j);
M = size(p.C,1); 
XHTCACHE = struct(); 
A = p.C(:,:,j);
Np = eye(size(A)); Np(q,q) = 0; 

if p.PC.type == 3, % Gamma-matrix.
    SC_inv = diag(p.PC.S_inv(q,:,j)); % this is the diagonal for this row.   
    muC = p.PC.mu(q,:,j)';
else
    SC_inv = p.PC.S_inv;
    muC = p.PC.mu;
end

b = SC_inv * muC;    


S_inv = SC_inv;
% S_inv
if p.debug,
    s0 = -1/2 * (A(q,:)'-muC)' * SC_inv * (A(q,:)'-muC);        
    ep = 1*((1:M)' == q);
    P = ep * ep';
    N = eye(M) - P;
    K = muC' * SC_inv * muC;
end
for k=1:length(p.s),
    ps = p.s(k);
    if size(Cnj,3) > 0, 
        dm = mmx('mult',Cnj,ps.s(:,1:end-1) );
    else 
        dm = zeros(size(Cnj,1),size(ps.s,2)-1); 
    end
    xh = ps.s(:,2:end) - sum(bsxfun(@times, dm,  permute(ps.v(1:end ~= j, 2:end), [3, 2, 1])),3) ;
    xh = bsxfun(@minus, xh, p.Pw.mu ) -  p.U * p.s(k).v(:,2:end);
    xt = bsxfun(@times, ps.v(j,2:end), ps.s(:,1:end-1));    
    XHTCACHE(k).xh = xh;
    XHTCACHE(k).xt = xt;    
    S_inv = S_inv + p.Pw.S_inv(q,q) * (xt * xt');
    b = b + xt * (xh' - xt' * A' * Np) * p.Pw.S_inv(:,q); 
    if p.debug,
        for t=1:size(xh,2),
            s0 = s0 - 1/2 * (xh(:,t) - A * xt(:,t))' * p.Pw.S_inv * (xh(:,t) - A * xt(:,t));
        end
        K = K + trace(xt(:,:,1)' * (N*A)' * p.Pw.S_inv(:,:,1) * (N*A) * xt(:,:,1));            
        K = K + trace( xh(:,:,1)' * p.Pw.S_inv(:,:,1) * xh(:,:,1)) ...
            -2* trace(N * A * xt(:,:,1) *xh(:,:,1)' * p.Pw.S_inv(:,:,1) );
    end
end    
if p.debug,
    po = p;
    ap = A(q,:)';
    res2 = -1/2 * (ap - S_inv\b)' * S_inv * (ap - S_inv\b) + 1/2 * (b' * (S_inv\b)) - 1/2 * K;
    if abs(res2-s0) > 10^-10,
        assert(false)
    end 
end 
A(q,:) = mvnrnd(S_inv\b,inv(S_inv));     
p.C(:,:,j) = A;  
 
for k=1:length(p.s),
    p.s(k).w(:,2:end) = bsxfun(@plus, p.Pw.mu, XHTCACHE(k).xh - A * XHTCACHE(k).xt);        
end    
if p.debug,
    MDS_selfcheck(p);
    MDS_dlp_test(p,po,A(q,:)',ap,S_inv\b,S_inv);
    %%
end
p.C(:,:,j) = A;   
p = MDS_fward(p);
end



