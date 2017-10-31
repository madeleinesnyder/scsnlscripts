%%
% Normal-half-cachy distribution. 
% Implemented for n observations from the normal distribution (this is not
% a multivariate normal distribution). 
% Implements the model
% x ~ N(0,sigma)
% sigma ~ HalfC(scale)

classdef th_NhalfC
    properties
        psigma
        sigma
        S_inv
        name
    end
    
    methods
        function p = th_NhalfC(scale)
            if nargin < 1, 
                scale = 10;  
            end            
            p.psigma = th_cauchy(scale);    
            p.sigma = abs(mvnrnd(0,p.psigma.S_inv));
            p.S_inv = 1/p.sigma^2; 
            x = mvnrnd(0,p.sigma );                    
            p.MCMC(x,1);            
            p.name = 'Normal-HalfCauchy';
        end
        
        
        function p = MCMC(p,D,samples)
            % mcmc sample for a while. 
            assert(size(D,2) == 1); 
            n = size(D,1);
            % p(x | p,a,b) = (a/b)^(p/2)/2/besselk(p,sqrt(a*b))*x^(p-1)*exp(-(a*x + b/x)/2)
            %%
            debug = false;
            if debug,
                close all;
                sigs = linspace(0.01, 10, 1000); 

                for i=1:length(sigs)
                    y(1,i) = limvnpdf(sigs(i), 0, p.psigma.S_inv);
                    y(1,i) = y(1,i) + limvnpdf(D,0, 1/sigs(i)^2);
                end

                y = exp(y-max(y)); y = y / trapz(sigs,y);
                plot(sigs,y)
            end
            %
            %p.psigma.S_inv                        
            %
            pp = -(n+1)/2 +1; 
            a = p.psigma.S_inv;
            b = sum(D.^2,1);
            
            z = gigrnd(pp,a,b,1);
            p.sigma = sqrt(z);
            %p.psigma.x = sigma;
            p.S_inv = 1/p.sigma^2; 
            p.psigma.MCMC(p.sigma); 
            
            
            
            if debug
                [v,xs] = hist(sigma,100);
                v = v/trapz(xs,v);
                hold all;
                plot(xs,v)
            end                        
            %p.psigma
            %sigma            
            %%                        
            %beta = sum(D.^2,1)/2 + p.psigma.beta
            %alpha = p.psigma.alpha            
            %v = lgamrnd(alpha,beta)            
            %%
            
        end
        
    end    
end
