function p = MDS_sample_hypers(p,opts)
if nargin < 1, test(); return; end
% HS = 0; 
% dbug = p.debug;
if p.debug, 
    lp(1) = MDS_logp(p);    
    p0 = p;
end 
[M,~,J] = size(p.C);
%%
for j=1:J,
    p.PC(j) = p.PC(j).MCMC(vec(p.C(:,:,j))' );
end
%%
if any(p.PC(1).sigma > 1000),
%     234
disp('bad, too large sigma!');
    assert(false);
end
% return
S = length(p.s);
[M,T] = size(p.s(1).e);

for m=1:M,
    E = zeros(S,T);
    W = zeros(S,T);
    SF = zeros(S,T);    
    for s=1:S,
        E(s,:) = p.s(s).e(m,:);
        W(s,:) = p.s(s).w(m,:);        
        SF(s,:) = p.s(s).sf(m,:);
    end
    p.Pe(m) = p.Pe(m).MCMC(E,4);
    if p.Pw_sigma_temporal,
        p.Pw(m) = p.Pw(m).MCMC(W,4);    
    end    
    %%
    if opts.Tsf > 0, 
        p.Psf(m) = p.Psf(m).MCMC(SF,4);
    end    
end

for m=1:M,
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
    %%
    if p.debug,
        s1 = -1/2 * mu0' * S_inv * mu0 + q' * mu0 -  1/2 * C;             
        if ~eqc(s0,s1)
           assert(false);
        end        
    end
    mus = imvnrnd(S_inv\q, S_inv);    
    if p.debug,
        p1 = p; 
        p1.Pe(m).mu(:) = mus; 
        MDS_dlp_test(p1,p,mus,p.Pe(m).mu(1),S_inv\q,S_inv);
        %%
        MDS_selfcheck(p);    
    end
    p.Pe(m).mu(:) = mus*1; 
end 

if p.debug,
    lp(end+1) = MDS_logp(p);    
end
3; 
%%
if ~p.Pw_sigma_temporal,
    Dw = [p.s.w]';
    p.Pw = p.Pw.MCMC(Dw);
end

%%
if p.debug,
    dlp3 = MDS_logp(p) - lp(end); 
end 
p = MDS_fward(p);
end