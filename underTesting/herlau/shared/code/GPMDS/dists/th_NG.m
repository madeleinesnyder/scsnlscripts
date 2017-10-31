% this is a prior of 
% sigma^-2 ~ Gam(a,b)
% x ~ N(0,sigma^2)     

classdef th_NG,
    
    properties 
        S_inv
        S
        L
        name
        gptype
        zero_lock
        alpha 
        beta
        mu
        s0
        tau
        sigma
        M
        type
    end
    methods
        
        function GP = tau2P(GP,tau2)
            GP.tau = tau2;
            GP.sigma = 1./sqrt(GP.tau);
            GP.S = GP.s0.^2 * GP.sigma.^2;
            GP.S_inv = diag(1./GP.S);
            GP.S = diag(GP.S);
            GP.L = cholcov(GP.S);

        end
        
        
        function P = th_NG(mu0,alpha,beta,s0) % how does the scaling work?
            if nargin < 1,
                M = 2;                 
                alpha = 1; 
                beta = 2; 
                s0 = 0.5; 
                mu0 = zeros(M,1);
                p = th_NG(mu0,alpha,beta,s0);
                 %%
                N = 1000; 
                sig = zeros(M,N);
                D = zeros(N,M);
                [~, sig] = p.sample();                
                for i=1:N,
                    D(i,:) = p.mu + randn(p.M,1) .* sqrt( diag(p.S) ); %p.sample()
                end
                hist(D(:,1), 20); hold on;
                hist(D(:,2), 20); hold on;
                T = 100; 
                ss = zeros(T,length(p.sigma));
                for t=1:T
                    [p,ss(t,:)] = p.MCMC(D,1);                
                end
                mean(ss,1)
                return;                 
            end  
            P.type = 100;
            P.name = 'Normal-Gamma (Fixed alpha-beta)';
            if nargin < 5, s0 = 1; end
            % with fixed alpha, beta?
            P.M = size(mu0,1);
%             P.M = M;
            tau2 = zeros(P.M,1) + alpha./beta;
            
%             P.sigma = 1./sqrt(alpha/beta + zeros(P.M,1)); 
%             tau2 = 1 ./ P.sigma.^2;                     
            P.s0 = s0; 
            P = P.tau2P(tau2); 
            
            P.mu = mu0;
            P.alpha = alpha;
            P.beta = beta;
            
            
        end
        
        
        function x = GP2vec(GP)
            x = [GP.l,GP.sigf,GP.sign];
        end        
        
        function GP = vec2GP(GP,x,doinvert)
            
            if nargin < 3, doinvert = true; end 
            GP.l = x(1); 
            GP.sigf = x(2);
            GP.sign = x(3);
             
            GP.S = K(GP, GP.x); % generic-ass kernel function
            if GP.zero_lock, mu = [0 ; GP.mu ; 0]; 
            else mu = GP.mu; end 
            GP.mu = mu; 
            if GP.zero_lock,
                [MU,SS] = GPmarg(GP, mu', [2:size(mu,1)-1 ]);
                GP.mu = MU'; 
                GP.S = SS; 
            end
            if doinvert,  
                GP.S_inv = inv(GP.S);                       
            else
                GP.S_inv = -inf;
            end
            if GP.sign == inf,
                GP.S_inv = zeros(size(GP.S_inv));
            end
            GP.L = cholcov(GP.S);            
        end
         
        function y = logp(GP,D)
            y = lmvnpdf_chol(D,GP.mu,GP.L);            
            y = y + sum(lgampdf( GP.tau, GP.alpha, GP.beta));
%             if GP.sigf > 0, y = y +  lgampdf(GP.beta); end
%             if GP.sign > 0, y = y + lgampdf(GP.alpha);  end
        end
        function [y,sigma] = sample(GP)            
            tau2 = gamrnd(GP.alpha,GP.beta,[GP.M,1]);
            sigma = 1./sqrt(tau2);
            y = normrnd(GP.mu,sigma);
        end
        
        function [GP,X] = MCMC(GP,D, nsamples)
            %% make basic mcmc RW sampler.
            if nargin < 3, nsamples = 1; end
            %%
            N = size(D,1);
            a2 = GP.alpha + N/2;
            b2 = GP.beta + 1/2 * sum(bsxfun(@minus, D, 0*mean(D,1)).^2,1) / GP.s0^2;             
            tau2 = gamrnd(a2, 1./b2);
            %%
            GP = GP.tau2P(tau2); 

            X = GP.sigma';
            
            %%
            return;
            %%
            
            x0 = log(GP.GP2vec());
            proprnd = @(x) x + randn(size(x))/10;
            I = isfinite(x0);
            
            function y = mylogp(lx)
                GP2 = GP; 
                dx = x0;
                dx(I) = lx;
                GP2 = GP2.vec2GP(exp(dx),false); 
                y = GP2.logp(D) + sum(lx(isfinite(lx))); 
            end
            X = mhsample_fast(x0(I),nsamples,'logpdf',@mylogp,'proprnd',proprnd,'symmetric',1);
            x0(end,I) = X(end,:);
            
            GP = GP.vec2GP(exp(x0));
            GP.S_inv = inv(GP.S);
        end
    end
end