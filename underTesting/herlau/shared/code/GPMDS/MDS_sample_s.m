function [p,S_inv_w,S_invY,q1,q2] = MDS_sample_s(p,only_matrices)
if nargin < 2, only_matrices = false;
if nargin < 1, test(); return; end

for sub=1:length(p.s),
    % construct observational matrix.
    %% construct the y-matrix.
    [M,T] = size(p.s(sub).y);
    if ~(isfield(p.s(sub),'c') && isfield(p.s(sub).c,'Spr')),  
       [M,T] = size(p.s(sub).y);
       pp = M; q = M; 
       r = 1; ss = T; 
       p.s(sub).c.Spr = PS(pp,r);
       p.s(sub).c.Sqs = PS(q,ss);  
    end     
    [~,q1,S_inv_w] = MDS_mk_S_matrix(p,sub,false);
    if false,        %% marginalize. 
        S_inv = S_inv_chol' * S_inv_chol;
        mu = S_inv \ q1
        t = 117; 
        J = t+1:length(mu);
        S = inv(S_inv);
        J1 = []; 
        J2 = []; 
        J3 = [];
        
        for m=1:M,
            J2 = [J2, (1:t-1) + (m-1)*T];
            J1 = [J1, t + (m-1)*T];
            J3 = [J3, (t+1:T) + (m-1)*T];
        end
        S11 = S(J1,J1);
        S12 = S(J1,J2); S21 = S12'; 
        S22 = S(J2,J2);
        mu1 = mu(J1,1);
        mu2 = mu(J2,1);
        
        mubar = mu1 + S12 * (S22\(mu2-mu2));
        Sbar = S11 - S12 * (S22\S21);        
        diag(Sbar)' -p.Pw.sigma2        
    end    
    
    [M,T] = size(p.s(sub).y);
    y = p.s(sub).y;
    H = p.s(sub).B * p.s(sub).Phi;

    S_invY = zeros(T*M,T*M);     
    L = size(H,2); 
    %%
    S_invY2 = zeros(T*M);
    q2 = zeros(T*M,1);
    for k=1:M,
        %%
        HR = zeros(T,T);
        for t=1:T,
            Lmin = min(t,L);
            HR(t, t-Lmin+1:t) = H(k,Lmin:-1:1); %end:-1:end-Lmin+1 ); %Lmin);        
        end
        ym = y(k,:)' - p.Pe(k).mu - p.s(sub).yhatf(k,:)';        
        I = T*(k-1)+1:T*k;
        
        S_invY2(I,I) = HR' * p.Pe(k).S_inv * HR;
        q2(I,:) = HR' * p.Pe(k).S_inv * ym;        
%       
    end
%    
    %%
    if false,
        %%
        clc
        C = p.C
        sigw = sqrt(p.Pw.sigma2)         
        S_inv_w(1:10,1:10)
        T = size(S_inv_w,1);
        L = zeros(T);
        for t=2:T,
            L(t,t-1) = 1; 
        end
        L(1:10,1:10); 
        dx  = 1/sigw * (-L *(1+C) + eye(T));
        sig = dx' * dx;
        sig(1:10,1:10)
        
        inv_sig = inv(sig); 
        diag(inv_sig(1:2:end,1:2:end))
        S = sigw.^2 / 2 * 1/abs(C) 
%         S = 1./S
        
%         1 / (C * sigw^2 )
    end
    %%
    S_inv = S_inv_w + S_invY2;
    qq = q1+q2;
    if rcond(S_inv) < 10^-20,
        disp('sample-s, bad rcond');%7888
    end
    muB = S_inv \ qq;        
    
    
    if ~only_matrices,
        %%
%         rcond(S_inv)
        if p.debug && (rcond(S_inv) < 10^-20 || isnan(rcond(S_inv))),
            assert(false);
        end
        s = imvnrnd(muB, S_inv)';    
        
        muB = reshape(muB,[T,M])';
        s = reshape(s,[T,M])';
        if isfield(p.s, 's_true'),
            HH = M+1; WW = 1; 
            for k=1:M,
                subplot(HH,WW,k);hold off;
                plot(muB(k,:),'k-');hold on;
                plot([s(k,:) ; p.s(sub).s_true(k,:); p.s(sub).y(k,:) ]');             
            end
        end
        if p.debug,
            p0 = p;
            p1 = p;
            p1.s(sub).s = s;
            p0 = MDS_fward(p0);
            p1 = MDS_fward(p1);
            MDS_dlp_test(p1,p0,vec(p1.s(sub).s'), vec(p0.s(sub).s') ,S_inv \ (q1+q2),S_inv);            
        end
        p.s(sub).s = s;     
    end
end
if ~only_matrices,
    p = MDS_fward(p);
end
end