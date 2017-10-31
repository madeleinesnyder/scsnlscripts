%% Simulates 
% x - C(s)
% using the Normal-Gamma represenatation of x where.
%
% x   ~ N(0,sig^-2=tau)
% tau ~ G(1/2, s^2/2) 
classdef th_cauchy < th_NG
    properties
    %x % the current value of the cauchy distribution. 
    end
methods 
    
    function p = th_cauchy(scale)
        if nargin < 1            
            scale = 10;
            234;
        end            
        alpha = 1/2; 
        beta = scale^2/2; 
        
        p = p@th_NG(0,alpha,beta);        
        x = 1;          
        p = p.MCMC(x);        
        p.name = 'Cauchy distribution'; 
    end    
    function p = MCMC(p,D)
        % MCMC sample the internal representation of x. 
        assert(size(D,2) == 1);        
        p = MCMC@th_NG(p,D);         
    end    
end

end
