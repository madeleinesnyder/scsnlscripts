classdef th_GAM,
% Implements a Gamma distribution. Idea is that $x ~ Gam(alpha,beta) and alpha,beta controlled
% by their own Gamma distributions.  During sampling (MCMC), alpha,beta are
% sampled. given $N$ observations of $x$ which must be $N \times 1$ (02450
% format). 
    properties
        alpha
        beta         
        Palpha
        Pbeta
    end
    
    methods
        function p = th_GAM(alpha,beta,opts)
            
            if nargin < 3, opts = struct(); end
            if nargin < 2, alpha = 1; end
            if nargin < 1, beta = 1; end
                
            p.alpha = alpha;
            p.beta = beta; 
            
%             p.Palpha.alpha = 1; 
%             p.Palpha.alpha = 1;
%             p.Pbeta.beta = 1;
%             p.Pbeta.alpha = 1;
            
            od.name = 'Gamma distribution (alpha,beta / wiki parameterization)'; 
            od.Palpha.alpha = 1;
            od.Palpha.beta = 1;
            od.Pbeta.alpha = 1;
            od.Pbeta.beta = 1;
            
            od.Palpha.sample = false;
            od.Pbeta.sample = false;
            
            opts = ssfr(od,opts);
            
            p.Palpha = opts.Palpha;
            p.Pbeta = opts.Pbeta;
                        
            if nargin < 1,
                close all;
                alpha0 = 1; beta0 = .10; 
                x = gamrnd(alpha0, 1/beta0, [100, 1]);
                % what exactly are we doing here?
                
                alpha0 = 1; beta0 = 1; 
                opts = struct();
                p = th_GAM(alpha0,beta0,opts);
                p.logp(x)    
                p.MCMC(x)
                
                234;
            end
        end
        function [lp,lpvec,lprest] = logp(p,x)
            lpvec = lgampdf(x,p.alpha,p.beta);
            lprest = 0;
            if p.Palpha.sample
                lprest = lprest + lgampdf(p.alpha, p.Palpha.alpha, p.Palpha.beta);
            end
            if p.Pbeta.sample
                lprest = lprest + lgampdf(p.beta, p.Pbeta.alpha, p.Pbeta.beta);
            end
            lp = sum(lpvec)+lprest;
        end
        
        function p = MCMC(p,X)
            %keyboard        
            if ~p.Palpha.sample && ~p.Pbeta.sample, return; end
            %%
            %close all;
            bb = linspace(0.01, mean(1./X)*2, 1000);
            for i=1:length(bb)
                break;
                p2 = p;
                p2.beta = bb(i);
                lp(i) = p2.logp(X);
            end
            %plot(bb,exp(lp-max(lp)));
            
            %%
            assert(size(X,2) == 1);            
            dd = .3; 
            sig_alpha = dd  * p.Palpha.sample;
            sig_beta = dd * p.Pbeta.sample;            
            
            %%
            function lp = mlp(vx)
                p2 = p;
                p2.alpha = exp(vx(1));
                p2.beta = exp(vx(2));                
                lp = p2.logp(X) + sum(vx);
            end
            function x2 = proposal(vx)                
                x2 = vx + randn(size(vx)) .* [sig_alpha, sig_beta];                 
            end
            t0 = tic();
            [s,accept] = mhsample( log([p.alpha,p.beta]), 1000, 'logpdf', @mlp,'proprnd',@proposal,'symmetric',1);
            t0 = toc(t0);
            fprintf('dt = %g, acc=%g\n', t0, accept);
            
%             plot(exp(s))
%             s(end,:)
            ss = exp(s(end,:)); 
            p.alpha = ss(1); p.beta = ss(2);     
            %keyboard;
        end        
    end    
end