function p = MDS_sample_U(p)
if nargin < 1, test(); return; end
for k=1:length(p.s),
    [~,T] = size(p.s(k).y);
    for h=1:length(p.PU), 
        if p.U_meanstim
            ddU = -p.dT * sum(bsxfun(@times, permute(p.s(k).v, [3, 2, 1]), mmx('mult', p.C, sum( bsxfun(@times, p.s(k).U(:,:,h), p.s(k).u(:,:,h)),3)) ), 3);
        else
            ddU = bsxfun(@times, p.s(k).U(:,:,h), p.s(k).u(:,:,h));
        end
        ds = p.s(k).w + ddU;
        %%
        if p.debug,
            p1 = p; 
            p1.s(k).U(:,:,h) = 0; 
            p1 = MDS_fward(p1);
            assert(eqc(ds, p1.s(k).w));            
        end
        %%
        M= size(p.C,1);
        
        if p.Pw_sigma_temporal,
            S_inv = p.PU(h).S_inv;
            q = zeros(M,1);
            for m=1:M, 
                if p.U_meanstim,
                    %%
                    dCmu = p.dT * sum( bsxfun(@times, bsxfun(@times, permute( p.C(m,:,:), [2, 1, 3]), p.s(k).u(:,:,h)), permute(p.s(k).v, [3, 2, 1]) ) , 3);                    
                    S_inv = S_inv + dCmu * p.Pw(m).S_inv * dCmu';
                    q = q - dCmu * p.Pw(m).S_inv * ds(m,:)';
                end
            end 
            if ~p.U_meanstim,
                %%
                SS = cell2mat(arrayfun( @(j)diag( p.Pw(j).S_inv),1:M,'UniformOutput',false))';                
                S_inv =  p.PU(h).S_inv + diag(sum(p.s(k).u(:,:,h).^2  .* SS,2));                
                q = sum((p.s(k).u(:,:,h).*ds) .* SS, 2);
            end
        else
            %%
            if p.U_meanstim,
%                 clc
                Cmu = zeros(M,T);
                S_inv = p.PU(h).S_inv;
                q = zeros(M,1);
                %
                L = cholcov(p.Pw.S_inv);
%                 p.Pw.S_inv - L' * L
                %
                for m=1:M,
                    dCmu = -p.dT * sum( bsxfun(@times, bsxfun(@times, permute( p.C(m,:,:), [2, 1, 3]), p.s(k).u(:,:,h)), permute(p.s(k).v, [3, 2, 1]) ) , 3);                    
                    Cmu(m,:) = p.s(k).U(:,:,h)' * dCmu;
                    Wt = kron(dCmu', L(:,m) ); 
                    S_inv = S_inv + Wt' * Wt;
                    q = q + vec( (ds' * p.Pw.S_inv(:,m))' *  (dCmu'));
                end
                %
                Cmu2 = zeros(M,T);
                for t=1:T,
                    %
                    Ct = sum(bsxfun(@times, p.C, permute(p.s(k).v(:,t),[2, 3, 1])), 3);
                    Cmu2(:,t) = -p.dT* Ct * p.s(k).u(:,t,h);
                end
                if p.debug,
                     s0 = -1/2* trace( diag(p.s(k).U(:,:,h)) * p.PU(h).S_inv * diag(p.s(k).U(:,:,h)) ) ...
                    -1/2 * trace(p.s(k).w' * p.Pw.S_inv * p.s(k).w);
                    
                     s1 = -1/2 * p.s(k).U(:,:,h)' * S_inv * p.s(k).U(:,:,h) + q' * p.s(k).U(:,:,h) - 1/2 * trace(ds' * p.Pw.S_inv * ds);
                     UU = p.s(k).U(:,:,h);
                     dsw = (ds - bsxfun(@times,  Cmu, UU) );
                     
                     
                     s2 = -1/2 * trace( diag(p.s(k).U(:,:,h)) * p.PU(h).S_inv * diag(p.s(k).U(:,:,h)) ) ...
                    -1/2 * trace( dsw' * p.Pw.S_inv * dsw);
                    
%                      [s0,s1,s2]-s0
                end
                345;
                %%
            else

                S_inv = p.PU(h).S_inv + diag(sum(p.Pw.S_inv * p.s(k).u(:,:,h) .* p.s(k).u(:,:,h),2));
                q = sum(bsxfun(@times,  p.Pw.S_inv * ds, p.s(k).u(:,:,h)), 2);
                
            end
            
            %%
            if false,
                clc
                s0 = -1/2* trace( diag(p.s(k).U(:,:,h)) * p.PU(h).S_inv * diag(p.s(k).U(:,:,h)) ) ...
                    -1/2 * trace(p.s(k).w' * p.Pw.S_inv * p.s(k).w);
                %                 s1 = trace( diag(p.s(k).U(:,:,h)) * p.PU(h).S_inv * diag(p.s(k).U(:,:,h)) );
                
                s1 = -1/2 * p.s(k).U(:,:,h)' * S_inv * p.s(k).U(:,:,h) + q' * p.s(k).U(:,:,h) - 1/2 * trace(ds' * p.Pw.S_inv * ds);
                
                dsw = (ds - diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h) );
                s2 = -1/2 * trace( diag(p.s(k).U(:,:,h)) * p.PU(h).S_inv * diag(p.s(k).U(:,:,h)) ) ...
                    -1/2 * trace( dsw' * p.Pw.S_inv * dsw);
                
                s3 = -1/2 * trace( diag(p.s(k).U(:,:,h)) * p.PU(h).S_inv * diag(p.s(k).U(:,:,h)) ) ...
                    -1/2 * trace( ds' * p.Pw.S_inv * ds) - 1/2 * trace( (diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h))' * p.Pw.S_inv * (diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h))) ...
                    + trace( ds' * p.Pw.S_inv * diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h));
                
                ss0 = -1/2 * trace( ds' * p.Pw.S_inv * ds);
                s4 = ss0 -1/2 * p.s(k).U(:,:,h)' * p.PU(h).S_inv * p.s(k).U(:,:,h) ...
                    - 1/2 * trace( (diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h))' * p.Pw.S_inv * (diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h))) ...
                    + trace( ds' * p.Pw.S_inv * diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h));
                
                - 1/2 * trace( (diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h))' * p.Pw.S_inv * (diag(p.s(k).U(:,:,h)) * p.s(k).u(:,:,h))) ...
                    + 1/2 * p.s(k).U(:,:,h)' * (S_inv-p.PU(h).S_inv) * p.s(k).U(:,:,h)
                
                [s0,s1, s2, s3, s4] - s0
            end
            
        end
        %             end
        U2 = imvnrnd(S_inv\q, S_inv)';
        
        if p.debug,
            %%
            p1 = p;
            %%
            p1.s(k).U(:,:,h) = U2;
            p1 = MDS_fward(p1);
            MDS_dlp_test(p1,p,U2,p.s(k).U(:,:,h),S_inv\q,S_inv);
        end

        p.s(k).U(:,1,h) = U2;
        p = MDS_fward(p);
        
    end
end

end
