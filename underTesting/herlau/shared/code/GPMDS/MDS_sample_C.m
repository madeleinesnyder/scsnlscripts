function p = MDS_sample_C(p)
if nargin < 1, test(); return; end
[M,~,J] = size(p.C);
for j=1:J,
    q_all = zeros(M*M,1); 
    S_inv_all = zeros(size(p.PC(j).S_inv));    
    
    D = S_inv_all; 
    B = D*0; 
    for sub=1:length(p.s),
       [M,T] = size(p.s(sub).s);
       %%
       if p.debug && false,
           dXA = zeros(M,T-1); 
           dXX = zeros(M,T-1);
           for t=1:T-1,
               dXA(:,t) = p.s(sub).s(:,t+1)-p.s(sub).s(:,t) - p.s(sub).U * p.s(sub).du(:,t+1);               
               for jj=1:J, 
                   if jj==j, continue; end
                   dXA(:,t) = dXA(:,t) - p.C(:,:,jj) * p.s(sub).s(:,t) * p.dT * p.s(sub).v(jj,t+1);
               end
               dXX(:,t) = p.s(sub).v(j,t+1) * p.s(sub).s(:,t) * p.dT;               
           end
       end         
       
        R = lagm(T);

        if p.U_meanstim,
            dUU = -p.dT * sum(bsxfun(@times, p.s(sub).U, p.s(sub).u),3);
        else
            dUU = sum(bsxfun(@times, p.s(sub).U, p.s(sub).u ),3 );
        end
        
        %%
        A = p.s(sub).s-p.s(sub).s*R; 
                        
        if p.U_meanstim,
            SRV = bsxfun(@minus, (p.s(sub).s * R ), sum(bsxfun(@times, p.s(sub).U, p.s(sub).u),3) ); 
        else
            SRV = (p.s(sub).s * R ); 
        end
        SRV = bsxfun(@times, SRV, permute(p.s(sub).v, [3, 2, 1])) * p.dT; %(jj,:) )
        for jj=1:J,
            if jj==j, continue; end
            A = A - p.C(:,:,jj) * SRV(:,:,jj);
        end
                
%         X = bsxfun(@times, (p.s(sub).s * R ), p.s(sub).v(j,:) )* p.dT;  
        X = SRV(:,:,j); 
        if p.U_meanstim,
%             X = X + dUU;
        else
            A = A - dUU; 
        end
        
        dd = A - p.C(:,:,j) * X;
%         dd - p.s(sub).w
        if p.debug,
           assert(eqc(dd,p.s(sub).w));
        end
        
        %%
       if p.Pw_sigma_temporal,
           %%
           Sdiag = cell2mat(arrayfun(@(Pw)diag(Pw.S_inv), p.Pw,'UniformOutput',false))';           
            
           S_inv = zeros(M*M);
           q = zeros(M*M,1);
           for m=1:M,
               em = 0 + ((1:M)' == m);  
               %% 
               dd = sqrt(diag(p.Pw(m).S_inv));
               XSd = bsxfun(@times, X, dd');
               %%
               S_inv = S_inv + kron( XSd * XSd', em * em' ); % this is same as prev. line but ensure symmetry.                
               q = q + vec((em * em') * A * p.Pw(m).S_inv * X'  );        
           end
           %%  
           if p.debug,
               Cv = vec(p.C(:,:,j));
                s0 = 0; 
                for kw=1:length(p.Pw),
                    s0 = s0 -1/2 * trace(dd(kw,:) * p.Pw(kw).S_inv * dd(kw,:)');
                end
                s1 = -1/2 * trace( (A .* Sdiag) * A') + q' * Cv - 1/2 * Cv' * S_inv * Cv;
%                 [s0, s1] - s0 
                 
                if ~eqc(s0,s1)
                    [s0, s1] - s0 
                    assert(false);
                end
            end           
       else
           S_inv = kron(X*X', p.Pw.S_inv);
           q = vec(p.Pw.S_inv *A*X');            

           B = B+S_inv;            
            if p.debug,
                Cv = vec(p.C(:,:,j)); %; Cv = Cv(:);
                s0 = -1/2 * trace(dd' * p.Pw.S_inv * dd);
                s1 = -1/2 * trace(A' * p.Pw.S_inv * A) + q' * Cv - 1/2 * Cv' * S_inv * Cv;
                [s0,s1]-s0;
                %%
                assert(eqc(s0,s1));
                
            end
        
       end
       
       q_all = q_all + q;       
       S_inv_all = S_inv_all + S_inv;        
    end 
    B= S_inv_all; 
    D = p.PC(j).S_inv;
    S_inv_all = S_inv_all + p.PC(j).S_inv; 
    rcn = rcond(S_inv_all);
    if rcn < 10^-20,  
%         assert(false);
%%
        B= S_inv_all; 
        D = p.PC(j).S_inv;
        L = cholcov(B);
        Dinv = diag(exp(-log(diag(D))));
        v = L * Dinv;
        n = length(L);
%         m = cholcov((eye(n) + v * D * v'));
%         [size(m), size(v)]
        %
        JDX = diag(D) > 10^4;
        Dm = diag(D);  Dm(JDX) = 1;  Dm = diag(Dm);
        Dp = diag(D)*0; Dp(JDX) = diag(D(JDX,JDX))-1; Dp = diag(Dp);
        %
        if rcond(B+Dm) < 10^-30,
           assert(false); 
        end
        L = cholcov(B + Dm);        
         
        for i=find(JDX)',
            x = sqrt(Dp(:,i));

            L = cholupdate(L,x,'+');
        end   
        
%         eqc(L' *L, D+B)
        if p.debug,
            %%
%             clc
            L0 = cholcov(B + D);
            L' * L -L0' * L0;
            L - L0;
            %%
            if ~eqc(L0,L), 
                assert(false)
            end 
        end       
        %%        
        
        C2vec = randn(1,size(L,1))/(L') +( L\(L' \ q_all) )';
        %%
    else 
        if rcond(S_inv_all) < 10^-19, 
            S_inv_all = S_inv_all + eye(size(S_inv_all,1))*10^(-12);  % numerical stability. 
        end
        if isnan(rcn) || rcn < 10^-20,
            assert(false);
        end 
        C2vec = p.PC(j).post_sample(vec(p.C(:,:,j)),S_inv_all,q_all);        
    end
    C2 = reshape(C2vec,[M,M]); 
   if p.debug, 
        p0 = p;  
        p1 = p;
        p1.C(:,:,j) = C2;
        p1 = MDS_fward(p1);     
        MDS_selfcheck(p);
        MDS_selfcheck(p1); 
        %%
%         clc
        %%
        p0 = MDS_fward(p0);     
        %%
        MDS_dlp_test(p1,p0,C2vec(:), vec(p0.C(:,:,j) ) ,S_inv_all \ q_all,S_inv_all);  
        %%
   end
   p.C(:,:,j) = C2;    
   p = MDS_fward(p);
end

end