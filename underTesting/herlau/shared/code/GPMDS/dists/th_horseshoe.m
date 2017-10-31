classdef th_horseshoe
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
        x
        s0
        tau
        sigma
        M
        type
        
        lambda
        
        gamma_lambda
        gamma_tau    
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
                
        function p = th_horseshoe(mu0) % how does the scaling work?
            if nargin < 1,
                M = 10; 
                mu0 = zeros(M,1);                 
                addpath('../');
                
                
                p = th_horseshoe(mu0);
                 %%
                N = 13; 
                %%
                p.S_inv 
                for i=1:N,
                    D(i,:) = imvnrnd(p.mu, p.S_inv);
                end
                nsamples = 100; 
                p.MCMC(D,nsamples);
                
                
                %%
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
            p.type = 101;
            % https://arxiv.org/pdf/1502.00560v2.pdf
            p.name = 'Horshoe shrinkage prior';
            if nargin < 2, t0 = 1; end
%             p.M = size(mu0,1);
            p.mu = mu0;
            p.gamma_lambda = 1; 
            p.gamma_tau = 1; 
            %%   
            if false,
                x = linspace(0.01,15);
                for i=1:length(x),
                    y(i) = lhalfcachypdf(x(i),0,p.gamma_lambda);            
                end            
                T = 10000; 
                for j=1:T,
                    yy(j) = halfcachyrnd(0,p.gamma_lambda);                
                end
                xx = linspace(min(x),max(x),200);
                [nn,ee] = histcounts(yy,xx);
                close all;
                plot(x,exp(y)/trapz(x,exp(y)),'k-'); hold all;
                xx2 = (ee(2:end)+ee(1:end-1) )/2;
                plot( xx2,nn/trapz(xx2, nn ),'ro-');            
            end
            %%
            p = p.sample();

        end
        function p = sample(p)
            p.lambda = halfcachyrnd(p.mu, p.gamma_lambda * ones(size(p.mu)));            
            p.tau = halfcachyrnd(0,p.gamma_tau);            
            p.S_inv = diag( (p.lambda.^2 * p.tau.^2).^(-1) );
                        
        end
        
        function x = GP2vec(GP)
            x = [GP.l,GP.sigf,GP.sign];
        end        
        
        function GP = vec2GP(GP,x,doinvert)
            
            if nargin < 3, doinvert = true; end 
            GP.l = x(1); 
            GP.sigf = x(2);
            GP.sign = x(3);
             
            GP.S = K(GP, GP.x); % generic kernel function
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
            
        end
          function s = Cplus_sample(p,X,gamma),
            alpha = 1/2; 
            beta = gamma^2/2;
            %N = length(X);
            [N,M] = size(X);
            
            a = alpha + M/2; 
            b = beta  + X.^2/2;

            s = gamrnd(a,1./b);
        end
        function [p,X] = MCMC(p,D, nsamples)
            %% make basic mcmc RW sampler.
            if nargin < 3, nsamples = 1; end
            %% how will this play out?
            % now do some sampling.
            %%
            ts = zeros(nsamples,1);
            ls = zeros(nsamples,length(p.lambda)); 
            function y = d1smp(ll, x0, nsamples)
                mylogp = @(x)ll(exp(x)) + x; 
                proprnd = @(x) x + randn()/10;
                X = mhsample_fast(log(x0), nsamples,'logpdf',mylogp,'proprnd',proprnd,'symmetric',1);
                y = exp(X(:,:));
            end
            for j=1:nsamples,
%                 tau_s = p.Cplus_sample(p.tau,p.gamma_tau);
%                 lambda_s = p.Cplus_sample(p.lambda,p.gamma_lambda)            ;

                [N,M] = size(D);
                
                for k=1:size(D,2),
                    % make likelihood function 
                    ll = @(lamj)sum(lnormpdf( D(:,k), 0, p.tau * lamj )) + lhalfcachypdf(lamj, 0, p.gamma_lambda);
                    lams = d1smp(ll,p.lambda(k),10);                
                   p.lambda(k) = lams(end);                
                end
                %
                ll = @(tau)sum(sum(lnormpdf( D', 0, p.tau .* p.lambda ))) + lhalfcachypdf(tau, 0, p.gamma_tau);           
                taus = d1smp(ll,p.tau,10);
                p.tau = taus(end);
            

                ts(j) = p.tau;
                ls(j,:) = p.lambda; 
            end
            %%
            close all;
            subplot(3,1,1);
            plot(ts); 
            subplot(3,1,2);
            plot(ls);             
            subplot(3,1,3);
            y = bsxfun(@times, ls, ts)
            plot(y)
            %%
            %%
            N = size(D,1);
            a2 = p.alpha + N/2;
            b2 = p.beta + 1/2 * sum(bsxfun(@minus, D, 0*mean(D,1)).^2,1) / p.s0^2;             
            tau2 = gamrnd(a2, 1./b2);
            %%
            p = p.tau2P(tau2); 
%             GP.tau = tau2; 
%             GP.sigma = 1./sqrt(GP.tau);
%             GP.S = GP.s0.^2 * GP.sigma.^2;
%             GP.S_inv = 1./GP.S;
%             GP.S = diag(GP.S);
            
            X = p.sigma';
            
            %%
            return;
            %%
            
            x0 = log(GP.GP2vec());
            proprnd = @(x) x + randn(size(x))/10;
            I = isfinite(x0);
            
            function y = mylogp(lx)
                GP2 = p; 
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