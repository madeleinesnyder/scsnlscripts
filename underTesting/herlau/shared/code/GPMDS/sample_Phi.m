function p = sample_Phi(p)
if nargin < 1, test(); return; end
for k=1:length(p.s),
    for pp=1:size(p.s(k).Phi,1),        
        p.s(k) = dsample_Phi(p,k,pp);
        if p.debug, MDS_selfcheck(p); end
    end    
end 
  
end
function ps = dsample_Phi(p,k,pp)
ps = p.s(k);
%% 
[M,T] = size(ps.y);
[np,L] = size(ps.Phi);

ep = 1*((1:np)' == pp);
N = eye(np) - ep * ep';
% p.PPhi.mu = p.PPhi.mu*0;
% p.PPhi.S_inv = p.PPhi.S_inv;

% 
yt = ps.y - permute( sum(bsxfun(@times, (ps.B)', mmx('mult',N * ps.Phi,ps.x)),1),[2,3,1]  ); 
xh = bsxfun(@times, ps.B(:,pp)', ps.x);

S_inv = p.PPhi.S_inv;
q = p.PPhi.S_inv * p.PPhi.mu ;

if p.debug, 
    C = p.PPhi.mu' * p.PPhi.S_inv * p.PPhi.mu;     
end

for m=1:M,
    ytm = bsxfun(@minus, yt(m,:), p.Pe(m).mu');
    xhm = permute( xh(:,m,:),[1,3,2] );    
    S_inv = S_inv + xhm * p.Pe(m).S_inv * xhm';
    q = q + xhm *  p.Pe(m).S_inv * ytm';    
    if p.debug,
        C = C + ytm * p.Pe(m).S_inv * ytm';
    end
end

if p.debug,
    Phip = ps.Phi(pp,:);
    s0 = -1/2 * (Phip' - p.PPhi.mu)' * p.PPhi.S_inv *  (Phip' - p.PPhi.mu);
    for m=1:M,
        s0 = s0 -1/2 * ps.e(m,:) * p.Pe(m).S_inv * ps.e(m,:)';
    end
    s1 = -1/2 * Phip * S_inv * Phip' + q' * Phip' - 1/2 * C;    
    if ~eqc(s0,s1),
        assert(false);
    end
end

if false,
%%
Sc_inv = p.PPhi.S_inv;
q = p.PPhi.S_inv * p.PPhi.mu;
for t=1:T,
    Sc_inv = Sc_inv + xh(:,:,t) * p.Pe.S_inv * xh(:,:,t)';
    q = q + xh(:,:,t) * p.Pe.S_inv * yt(:,t);
end

if p.debug,
    po = p; 
    phip = ps.Phi(pp,:)';    
    
    K = -1/2 * trace(yt' * p.Pe.S_inv * yt) - 1/2 * p.PPhi.mu' * p.PPhi.S_inv * p.PPhi.mu;
    s2 = -1/2 * phip' * Sc_inv * phip + q' * phip + K;
    s3 = -1/2 * (phip - Sc_inv \ q)' * Sc_inv * (phip - Sc_inv \ q) + 1/2 * q' * (Sc_inv\q) + K;    
    s0 = -1/2 * (phip - p.PPhi.mu)' * p.PPhi.S_inv * (phip - p.PPhi.mu);
    for t=1:T,
        s0 = s0 + -1/2 * (ps.y(:,t) - ps.yhat(:,t) - p.Pe.mu )' * p.Pe.S_inv * (ps.y(:,t) - ps.yhat(:,t) - p.Pe.mu );        
    end
    er = abs(s0-s3) / abs(s0+s3);
    if er > 10^-10,
        assert(false);
    end
end
%% sample:
end
ps.Phi(pp,:) = imvnrnd(S_inv \ q, S_inv);


if false,
    I = inv(Sc_inv);
    I = (I + I')/2;
    Sc_inv \ q;
    ps.Phi(pp,:) = mvnrnd(Sc_inv \ q, I );
end
%% now we got to update: yhat, e.
ps.yhat = permute( sum(bsxfun(@times, ps.B', mmx('mult',ps.Phi,ps.x)),1),[2,3,1]  );
ps.e = ps.y - ps.yhat;

if p.debug,
    p2 = p;
    p2.s(k) = ps;
    
    MDS_dlp_test(p2,p,ps.Phi(pp,:)',Phip',S_inv\q,S_inv);    
end 

end