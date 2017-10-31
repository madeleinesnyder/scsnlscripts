function p = MDS_sample_C_WORKING(p)
if nargin < 1, test(); return; end
[M,~,J] = size(p.C);
for j=1:J,
    q_all = zeros(M*M,1); 
    JJ = M*M * (j-1) + 1:M*M*j;
    S_inv_all = p.PC.S_inv(JJ,JJ);    
    
    D = S_inv_all; 
    B = D*0;
    for sub=1:length(p.s),
       [M,T] = size(p.s(sub).s);
       %%
       if p.debug || true,
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
       X = dXX;
       A = dXA;
       
%        SD_inv = p.Pw.S_inv; %diag(1./diag(p.Pw.D));
       dd = A - p.C(:,:,j) * X;   
       if p.Pw_sigma_temporal,
           s0 = -1/2 * trace(dd * p.Pw.S_inv * dd');
           
           
           345;
       else
           S_inv = kron(p.Pw.S_inv, X*X');
           q = vec(X * A' * p.Pw.S_inv);       
           q_all = q_all + q;       
           B = B+S_inv;       
           S_inv_all = S_inv_all + S_inv;
       end
       
       if p.debug,
           Cv = vec(p.C(:,:,j)'); %; Cv = Cv(:);
           s0 = -1/2 * trace(dd' * p.Pw.S_inv * dd);
           s1 = -1/2 * trace(A' * p.Pw.S_inv * A) + q' * Cv - 1/2 * Cv' * S_inv * Cv;
           eqc(s0,s1);
       end
    end 
    rcn = rcond(S_inv_all);
    if rcn < 10^-20 || false, 
        assert(false);
        L = cholcov(B);
        Dinv = diag(exp(-log(diag(D))));
        v = L * Dinv;
        n = length(L);
        m = cholcov((eye(n) + v * D * v'));
        [size(m), size(v)]
        %%
        JDX = diag(D) > 10^4;
        Dm = diag(D);  Dm(JDX) = 1;  Dm = diag(Dm);
        Dp = diag(D)*0; Dp(JDX) = diag(D(JDX,JDX))-1; Dp = diag(Dp);
        %%
        L = chol(B + Dm);        
         
        for i=find(JDX)',
            x = sqrt(Dp(:,i));

            L = cholupdate(L,x,'+');
        end                
        if p.debug,
            %%
            clc
            L0 = cholcov(B + D);
            L' * L -L0' * L0
            L - L0
            %%
            if ~eqc(L0,L),
                45
            end 
        end      
        %%
        C2vec = randn(1,size(L,1))/(L') +( L\(L' \ q_all) )';
        %%
    else 
        C2vec = p.PC.post_sample(vec(p.C(:,:,j)'),S_inv_all,q_all);        
    end

    C2 = reshape(C2vec,[M,M])'; 
   if p.debug, 
        p0 = p;  
        p1 = p;
        p1.C(:,:,j) = C2';% = s ; 
        p1 = MDS_fward(p1);     
        MDS_selfcheck(p);
        MDS_selfcheck(p1);
        
        MDS_dlp_test(p1,p0,C2vec, vec(p0.C(:,:,j)') ,S_inv_all \ q_all,S_inv_all);  
        %%
   end
   p.C(:,:,j) = C2;    
   p = MDS_fward(p);
end

end