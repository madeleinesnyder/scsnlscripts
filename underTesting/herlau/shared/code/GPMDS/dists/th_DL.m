classdef th_DL
    % Dirichlet-Laplace shrinkage prior. 
% see http://www.tandfonline.com.globalproxy.cvt.dk/doi/pdf/10.1080/01621459.2014.960967?needAccess=true    
 properties 
     type
     name
        S_inv
        mu
        s
        theta
        tau        
        a
        
        stats
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
                
        function p = th_DL(mu0,a)
            % only allow N=1 so far. Must think more closely about the
            % general case. 
            if nargin < 1,
                M = 10; 
                mu0 = zeros(M,1);                 
                addpath('../');
                a = 1/M; 
                p = th_DL(mu0,a);
                [x,p] = p.sample();
                %% run some MCMCM.    
                T = 1000; 
                sigma = zeros(T,M);                 
                theta = zeros(T,M);
                tau = zeros(T,M);
                for t=1:T,
                    [p,sigma(t,:)] = p.MCMC(x);                       
                    theta(t,:) = p.theta;
                    tau(t) = p.tau;
                    t
                end
                %%
                close all;
                subplot(2,1,1);
                plot( abs(x),'ko');
                hold on;
                plot(1:M, sigma','r.');                
                subplot(2,1,2);
                %plot(x, mean(sigma,1),'o');
                
                
                plot(abs(x),mean(sigma,1),'o');
                return;               
            end  
            p.type = 102;
            % https://arxiv.org/pdf/1502.00560v2.pdf
            p.name = 'DL_a prior';
            p.a = a;
            p.mu = mu0;
            %%
            [x,p] = p.sample();  
            p.stats.sigma = 1; 
        end
        function [x,p] = sample(p)
            M = size(p.mu,1); 
            p.tau = gamrnd(M*p.a,1 / (1/2));
            p.theta = dirrnd(zeros(1,M) + p.a,1);
            p.s = exprnd(1/(1/2), [1,M]); % matlabs stupid parameterization strikes again
            sig2 = p.s.^2 .* p.theta.^2 * p.tau;
            p.S_inv = diag(1./sig2);
            x = imvnrnd(p.mu,p.S_inv);            
        end  
        function y = lDE(p,x,tau)
            y = -log(2*tau) - abs(x)./tau;
        end 
        
        function y = logp(p,x)
           % this could be prettier...
           M = size(p.mu,1);
           y = limvnpdf(x,p.mu,p.S_inv);
           
           y = y + sum(exppdf(p.s, 2));
           y = y + dirpdf(p.theta,zeros(M,1)' + p.a);
           y = y + lgampdf(p.tau,M*p.a,1/2);
           
            
        end 
       function X2 = post_sample(p,X,S_inv_all,q_all)
            X2 = imvnrnd(S_inv_all \ q_all, S_inv_all);              
       end
        function [p,X] = MCMC(p,x)
            %%
            M = size(x,2);
            T = zeros(1,M);
            t0 = tic();
            for k=1:M, 
                T(k) = gigrnd(p.a-1,1,2*abs(x(k)-p.mu(k)),1);
            end  
            if any(T == 0)% == 0, 
                %%
%                 clc
                for k=find(T==0),
                    abs(x(k)-p.mu(k))
                    z0 = gigrnd(p.a-1,1,2*abs(x(k)-p.mu(k)),1)
                    bb = (2*abs(x(k)-p.mu(k)));
                    aa = 1;
                    z1 = gigrnd(p.a-1,aa*bb,1,1);
                    z1  = z1*bb
                    T(k) = z1; 
%                     keyboard;
                end
                
                
                
%                 324
            end 
            if toc(t0) > 3, 
                disp('BAAAD');
            end
            
            %%
            p.theta = T / sum(T);          
%             log(abs(x-p.mu'))
%             if any(abs(x-p.mu') < 10^-12),
%                 342
%             end 
            t0 = tic(); 
%             Q = sum( abs(x-p.mu')./p.theta ); 
%             log(p.theta)
%             fprintf('sampling tau... %g', Q);
            
            p.tau = gigrnd(M*p.a-M, 1, 2 *sum( abs(x-p.mu')./p.theta ), 1);
%             fprintf(' done!\n');
            
            if toc(t0) > 3, 
                disp('BAAAD 2');
            end
            nu = p.theta * p.tau ./ abs( x - p.mu');            
            shat = zeros(1,M);
            
            for m=1:M,
                %iG = makedist('InverseGaussian','mu',nu(m),'lambda',1);
                %shat(m) = iG.random();     
                shat(m) = igrnd(nu(m),1);            
            end
            p.s = 1./shat;            
            p.S_inv = diag( 1./(p.s .* (p.tau.*p.theta).^2) );            
            X = 1./sqrt(diag(p.S_inv));
            
        end
        
      
    end
    
end