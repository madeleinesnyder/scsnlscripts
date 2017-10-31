%%
% This assumes that W ~ p(W | S_inv, mu)
% where S_inv is diagonal is the prior for the model.
% and that W is restricted to lie in a region defined by the variable
% I.
%
% SEE:
% http://mlg.eng.cam.ac.uk/pub/pdf/Sch09b.pdf
classdef th_MVN_constrained
    properties
        mu
        sigma
        M
        IV
        opts
        name = 'th-mvn-constrained';
    end
    methods
        function pWc = th_MVN_constrained(mu,sigma,opts)
            if nargin < 1,
                K = 4;
                mu = zeros(K,1);
                sigma = ones(1,K)/2;
                
                IV = [0, 1 ; -inf,inf ; -.5, .5 ;  .4, .4];
                pWc = th_MVN_constrained(mu,sigma,IV);
                %%
                X2 = pWc.post_sample(mu', diag(1./pWc.sigma),mu)
                lp = pWc.log_p_W_given_Sigma(X2,sqrt(pWc.sigma));
                return;
            end
            pWc.mu = mu;
            pWc.sigma = sigma;
            %         if size(IV,1) > 0,
            %         pWc.eqC = IV(IV(:,1) == IV(:,2),:);
            %         pWc.ineqC = IV(IV(:,1) ~= IV(:,2),:);
            if ~isfield(opts,'IV'),
                opts.IV = zeros(0,2);
            end
            
            pWc.IV = opts.IV;
            pWc.opts = opts;
            if isfield(opts,'name'),
                pWc.name = opts.name;
            end
            
            %         end
            pWc.M = length(mu);
        end
        function lp = log_p_W_given_Sigma(pWc,X,sigma)
            %%
            if size(pWc.IV,1) == 0,
                I = pWc.IV(:,1) ~= pWc.IV(:,2);
            else
                I = true(length(pWc.mu),1);
            end
            lp = sum(lnormpdf(X(:,I), pWc.mu(I,1), sigma(1,I)') );
            %        [Qx,qx,Rx,rx] = p.intervals2Qq(p.IV);
        end
        function [Qx,qx,Rx,rx] = intervals2Qq(pWc,IV)
            %%
            qx = [];
            Qx = [];
            for j=1:size(IV,1),
                qx(end+1,1) = -IV(j,1);
                qx(end+1,1) = IV(j,2);
                
                Qx(1+(j-1)*2,j) = -1;
                Qx(2+(j-1)*2,j) = 1;
                if IV(j,1) == IV(j,2),
                    qx(end-1:end) = NaN;
                    Qx(end-1:end,:) = NaN;
                end
                
            end
            I = ~isnan(qx);%IV(:,1) ~= IV(:,2);
            
            Qx = Qx(I,:); qx = qx(I,:);
            Qx = Qx';
            
            %%
            M = size(IV,1);
            Rx = zeros(0,M);
            rx = zeros(0,1);
            for j=1:size(IV,1),
                if IV(j,1) == IV(j,2),
                    Rx(end+1,j) = 1;
                    rx(end+1,1) = IV(j,1);
                end
            end
            Rx = Rx';
            
        end
        
        % compute the marginal of x_I given x_{~I} = a. I is a binary vector.
        function [Sbar_inv_I, mubar_I] = marginal_distribution(p, S_inv, mu, a, I)
            if nargin < 1,
                %234
                M = 5;
                S_inv = rand(M,M); S_inv = S_inv' * S_inv;
                mu = rand(M,1);
                
                I = rand(M,1) < 0.5;
                nnz(I)
                a = 0*rand(nnz(~I),1);
                %%
            end
            if nnz(~I) == 0, 
                Sbar_inv_I = S_inv;
                mubar_I = mu;
                return;
            end
            Sbar_inv_I = S_inv(I,I);
            mubar_I = mu(I,:) - S_inv(I,I) \ (S_inv(I,~I)*(a-mu(~I,:)) );
            
            if false
                S = inv(S_inv);
                S11 = S(I,I);
                S22 = S(~I,~I);
                S12 = S(I,~I); S21 = S12';
                
                mubarB = mu(I,:) + S12 * inv(S22) * (a - mu(~I,:));
                SbarB = S11 - S12 * inv(S22) * S21;
                Sbar_inv_B = inv(SbarB);
            end
            %%
        end
        
        % X must be in format (M x 1)
        function X2 = post_sample(p,X,S_inv,q_all)
            X0 = X;
            %%
            % alternative idea: First get rid of the equality constraints. 
            % by exact marginalization & using they apply to 
            % axis-aligned coordinates. 
           
            %%
            if size(p.IV,1) > 0,
                %%
                X = X0;
                [Qx,qx,~,~] = p.intervals2Qq(p.IV);
                %%
                %                 mux = S_inv \ q_all;
                % this part will deal with the exact marginalization.
                % Notice the constraints are re-defined. 
                
              Imarg = p.IV(:,1) ~= p.IV(:,2); % these are to be marginalized out. 
              a = X(~Imarg,:); 
                
              [Sbar_inv_I, mubar_I] = p.marginal_distribution(S_inv, S_inv\q_all, a, Imarg);
                
              muy = mubar_I;
              Sy_inv = Sbar_inv_I;
              
                %
                if false,
                    [U,S,~] = svd(Rx);
                    %         T = U';
                    %%
                    if size(S,2) ~= size(S,1),
                        S(size(S,1), size(S,1)) = 0;
                    end
                    %%
                    IT = diag(S) ~= 0;

                    T = U(:,IT)';
                    Tort = U(:,~IT)';
                    x0 = (Rx')\rx;
                    M = size(Tort,2);
                end                
                y = X(Imarg,:);                     
                Qy = Qx(Imarg,:);                    
                qy = qx; 
                %% 
                L = cholcov(Sy_inv);
                if size(L, 1) <= 1
                    keyboard; 
                end
                z  = L * (y-muy);
                Qz = (L') \ Qy;
                qz = qy - Qy' * muy; 
                
                for i=1:size(z,1)
                    d = Qz(i,:)';
                    n = qz - Qz(1:end ~= i,:)' * (z(1:end ~= i,:) );
                    nd = (n ./ d);
                    l = max( [-inf ; nd(isfinite(nd) & d < 0) ] );
                    u = min( [ inf ; nd(isfinite(nd) & d > 0) ] );
                    
                    lu(i,:) = [l,u];
                    luinv = normcdf([l,u],0,1);
                    if luinv(1) ~= 1 && luinv(2) ~= 0
                        z(i) = norminv( rand() * (luinv(2)-luinv(1) ) + luinv(1) );
                    end
                    %% 
                    if luinv(1) > 1-10^-10
                        z(i) = l;                        
                    end
                    %%
                    if luinv(2) <10^-10
                        z(i) = u;
                    end
                    %%
                end
                %%
                X2 = L \ z + muy; 
                X2(~Imarg,:) = a; % = X; 
                X2(Imarg,:) =  L \ z + muy; 
                if max(X2(:)) > max(p.IV(:)),
                    234
                end
                %%
            else
                X2 = imvnrnd(S_inv \ q_all, S_inv);                
            end
        end
        
        function b = satisfy_constraints(p,X)
            if size(p.IV,1) == 0, b = true; return; end
            b1 = all(X(:) <= p.IV(:,2) & X(:) >= p.IV(:,1));
            [Qx,qx] = p.intervals2Qq(p.IV);
            b2 = all(Qx' * X(:) - qx <= 0) ;            
            assert(b1 == b2);
            b = b1; 
        end
    end
end

