function p = sample_hypers(p,opts)
if nargin < 1, test(); return; end
HS = 0; 
dbug = p.debug;
if dbug, 
    lp(1) = MDS_logp(p);    
    p0 = p;
end 
[M,~,J] = size(p.C);
De = [p.s.e]'; 
%%
p.PC = p.PC.MCMC(vec(permute(p.C, [2, 1, 3]))' );


 
S = length(p.s);
[M,T] = size(p.s(1).e);
for m=1:M,
    E = zeros(S,T);
    for s=1:S,
        E(s,:) = p.s(s).e(m,:);
    end
%     lp = lp + p.Pe(m).logp(E);
    p.Pe(m) = p.Pe(m).MCMC(E,4);
    
end


for m=1:M,
    %p.Pe(m) = p.Pe(m).MCMC(De',4);
    % this is not so hot.
    
    
    
    %%
    mu0 = p.Pe(m).mu(1);
    S_inv = 0;
    q = 0;
    C = 0;
    s0 = 0; 
    for k=1:length(p.s),
        dy = (p.s(k).y(m,:) - p.s(k).yhat(m,:))';
        S_inv = S_inv + sum(sum(p.Pe(m).S_inv)); 
        
        q = q + sum(dy' * p.Pe(m).S_inv);
        C = C + dy' * p.Pe(m).S_inv * dy;
        if p.debug,
            s0 = s0 -1/2 * (dy - mu0)' * p.Pe(m).S_inv * (dy-mu0);
        end
    end    
    if p.debug,
        s1 = -1/2 * mu0' * S_inv * mu0 + q' * mu0 -  1/2 * C;        
        if ~eqc(s0,s1),
            assert(false);
        end        
    end
    %%
     
%     MDS_selfcheck(p)
     
    mus = imvnrnd(S_inv\q, S_inv);
    p.Pe(m).mu(:) = mus;     
    if p.debug,
        MDS_selfcheck(p);
    
    end
end

if p.debug,
    lp(end+1) = MDS_logp(p);    
end
3; 
%%
Dw = [p.s.w]';
p.Pw = p.Pw.MCMC(Dw);

%%
if p.debug,
    dlp3 = MDS_logp(p) - lp(end); 
end 
p = MDS_fward(p);
end