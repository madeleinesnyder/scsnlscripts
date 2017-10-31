function [p,S_inv_w,S_invY,q1,q2] = MDS_sample_sf(p)
if nargin < 2, only_matrices = false; end
if nargin < 1, test(); return; end
   
for sub=1:length(p.s),
   [M,T] = size(p.s(sub).y);
    for m=1:M,
%        w34
    % construct observational matrix
    %% construct the y-matrix.
    [S_inv_chol] = cholcov(p.Psf(m).S_inv);
%     q1 = p.Psf(m).mu;
    q1 = (S_inv_chol' * S_inv_chol) * p.Psf(m).mu;
    
    [M,T] = size(p.s(sub).y);
    y = p.s(sub).y;
    H = p.s(sub).B * p.s(sub).Phi;

    S_invY = zeros(T,T);     
    L = size(H,2); 
    %%
    S_invY2 = zeros(T);
    q2 = zeros(T,1);

    k = m; 
    %%
    HR = zeros(T,T);
    for t=1:T,
        Lmin = min(t,L);
        HR(t, t-Lmin+1:t) = H(k,Lmin:-1:1); %end:-1:end-Lmin+1 ); %Lmin);        
    end
    ym = y(k,:)' - p.Pe(k).mu - p.s(sub).yhat(k,:)';        
    I = T*(1-1)+1:T*1;

    S_invY2(I,I) = HR' * p.Pe(k).S_inv * HR;
    q2(I,:) = HR' * p.Pe(k).S_inv * ym;        
       
    S_inv_w = S_inv_chol' * S_inv_chol;
    
    S_inv = S_inv_w + S_invY2;
    qq = q1+q2;
    muB = S_inv \ qq;        
    
    if false,        
        lp1 = -1/2 * trace(p.s(sub).w' * p.Pw.S_inv * p.s(sub).w ); 
        for k=1:M,
            lp1 = lp1 - 1/2 * (p.s(sub).e(k,:)-p.Pe(k).mu') * p.Pe(k).S_inv * (p.s(sub).e(k,:)-p.Pe(k).mu')'; % this is the real answer.
        end
        % now test factorization.
        Xv = vec(p.s(sub).s');
        UU = p.s(sub).U * p.s(sub).du;
        lp2 = -1/2 * Xv' * S_inv * Xv + qq' * Xv  -1/2 * trace( UU' * p.Pw.S_inv * UU);
        for k=1:M,
            lp2 = lp2 -1/2 * (y(k,:)-p.Pe(k).mu') * p.Pe(k).S_inv * (y(k,:)-p.Pe(k).mu')';
        end
        %%
        if ~eqc(lp1,lp2)
           assert(false);
        end  
    end    
    if ~only_matrices,
        sf = imvnrnd(muB, S_inv)';            
        sf = reshape(sf,[T,1])';
        if p.debug,
            p0 = p;
            p1 = p;
            p1.s(sub).sf(m,:) = sf;
            p0 = MDS_fward(p0);
            p1 = MDS_fward(p1);
            MDS_dlp_test(p1,p0,vec(p1.s(sub).sf(m,:)'), vec(p0.s(sub).sf(m,:)') ,S_inv \ (q1+q2),S_inv);            
        end
        p.s(sub).sf(m,:) = sf;     
    end
    end
end
if ~only_matrices,
    p = MDS_fward(p);
end
end