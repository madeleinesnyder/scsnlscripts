classdef th_GBOX
% Gaussbox
% implements 
% http://mlg.eng.cam.ac.uk/pub/pdf/Sch09b.pdf
   properties 
        Q
        q        
   end
    methods
        %% 
        function y = logp(BX,X),
            y = limvnpdf(X,p.mu,p.S_inv);
            y = y + sum(lgampdf(p.sigma,p.alpha,p.beta));     
            y = y - sum(3);
        end
        function [Q,q] = intervals2Qq(BX,I)
            q = [];
            Q = [];
            for j=1:size(I,1),
                q(end+1,1) = -I(j,1);
                q(end+1,1) = I(j,2);
                
                Q(1+(j-1)*2,j) = -1;
                Q(2+(j-1)*2,j) = 1;
                
            end
            Q = Q';
        end
        % sigma ~ Gam(alpha,beta)
        % variables restricted to mu0, axis-aligned box. 
        function BX = th_GBOX(mu0,alpha,beta,IV)
            % construct a GPBOX with the given intervals (as many as there
            % are dimensions) and with means, variances given by Sigma. 
            if nargin < 1,
                %%
                close all;
                M = 1; 
                Sig = eye(M)*.5 + (1-eye(M))*(-0.7);
                mu = zeros(M,1); 
                mu(1) = -.8;
                
                N = 100000; 
                Xs = mvnrnd(mu,Sig,N);
%                 plot(Xs(:,1), Xs(:,2),'k.');
                
                I = ones(1,1) * [-1,1];
                I(1,:) = [-2.5,0.01];
                if M > 1, 
                    I(2,:) = [-1, .5];
                end 
                
                alpha = 2; beta = 2; 
                mu0 = mu*0; 
                
                BX = th_GBOX(mu0,alpha,beta,I);
                
                I = sum( bsxfun(@minus, Q' * Xs', q) >= 0,1)' == 0;       
                if M == 2, 
                    subplot(2,2,1);
                    plot(Xs(:,1), Xs(:,2),'k.'); hold on;
                    plot(Xs(I,1), Xs(I,2),'r.');                         
                end
                
                x = zeros(M,1);
                T = nnz(I); 
                xs = zeros(M,T);
                for t=1:T,
                    [~,x] = BX.MCMC(mu,Sig,x);
                    xs(:,t) = x; 
                end
                xs = xs'; 
                Xs = Xs(I,:); %, Xs(I,2)
                if M == 2, 
                    subplot(2,2,2);                    
                    plot(Xs(:,1), Xs(:,2),'k.'); hold on;
                    subplot(2,2,3);
                    plot(xs(:,1), xs(:,2), 'k.'); 
                    4

                    %%
                    subplot(2,2,4);
                    [nn1,cc1] = hist3(Xs)
                    mg1 = meshgrid(cc1{1}, cc1{2})
                    I1 = nn1 / sum(nn1(:));

                    [nn2,cc2] = hist3(xs)
                    mg2 = meshgrid(cc2{1}, cc2{2})
                    I2 = nn2 / sum(nn2(:));

                    imagesc([I1,I2]); %nn / sum(nn) )
                else
                    %%
                    B = 20;
                    hold off;
                    [y1,x1] = hist(Xs,B);
                    [y2,x2] = hist(xs,B);
                    y1 = y1 / trapz(x1,y1);
                    y2 = y2 / trapz(x2,y2);
                    
                    plot(x1,y1,'k.-'); hold all;
                    plot(x2,y2,'r.-'); hold all;                    
                end
                return; 
            end
            
            [Q,q] = intervals2Qq(BX,I);
            BX.Q = Q;
            BX.q = q;
            
        end
        function [BX,xn] = MCMC(BX, mu, Sig, x),
            %%
            x = x(:);
            M = length(x);
            L = cholcov(Sig)'; 
%             L*L' - Sig
            z  = (L')\(x-mu);
            Qz = L * BX.Q;
            qz = BX.q - BX.Q' * mu;
            for i=1:M,
                d = Qz(i,:)';
                n = qz - Qz(1:end ~= i,:)' * (z(1:end ~= i,:) );
                nd = (n ./ d);
                l = max( [-inf ; nd(isfinite(nd) & d < 0) ] );
                u = min( [ inf ; nd(isfinite(nd) & d > 0) ] );
                
                luinv = normcdf([l,u],0,1);
                z(i) = norminv( rand() * (luinv(2)-luinv(1) ) + luinv(1) );                
            end
            
            xn = L' * z + mu;
        end
    end
end