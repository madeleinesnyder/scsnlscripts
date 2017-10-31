classdef th_GP,
    properties 
        S_inv
        S
        L 
        name
        sigf
        sign
        l
        gptype
        zero_lock
        mu
        x
        opts
        stats
    end
    methods
        % type == 1 : GP.
        % type == 2 : Ohrenstein-Uhlenbeck
        
        function GP = th_GP(gptype,mu0,l,sigf,sign,opts)
            if nargin < 6, 
                opts = struct();
            end
            od.zero_lock = false;
            if islogical(opts),
                od.zero_lock = opts;
                opts = struct();
            end
            if l == 0 || sigf == 0, 
                l = 0; sigf = 0; 
            end
            
            od.Tsigf = 1; 
            od.Tsign = 1; 
            od.Tl = 1; 
              
            od.Pl.alpha = 3/2; 
            od.Pl.beta = 3/2; 
            
            od.Psign.alpha = 1;
            od.Psign.beta = 2;
            
            od.Psigf.alpha = 1;
            od.Psigf.beta = 2;
            
            
            opts = ssfr(od,opts);
            
            GP.gptype = gptype;
            GP.name = '';
            GP.zero_lock = opts.zero_lock;
            T = numel(mu0);
            
            GP.mu = mu0; 
            GP.x = (1:T+2*GP.zero_lock)';
            GP.opts = opts;
            GP.stats.l = [];
            GP.stats.sign = [];
            GP.stats.sigf = [];
                        
            if gptype == 1, GP.name = 'Gaussian Process';
            else GP.name = 'Ohrenstein-Uhlenbeck'; end
            GP = vec2GP(GP,[l,sigf,sign]);                                            
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
        function fbar = posterior(GP,xstar,y,sign)
            % http://www.gaussianprocess.org/gpml/chapters/RW2.pdf
            % eq. 2.21.
            fbar = K(GP, xstar, GP.x) *( (K(GP,GP.x,GP.x) + eye(size(GP.mu,1)) * sign) \ y);            
        end 
        function y = logp(GP,D)
            y = lmvnpdf_chol(D,GP.mu,GP.L);    
            if GP.sigf > 0, y = y + lgampdf(GP.sigf,GP.opts.Psigf.alpha,GP.opts.Psigf.beta); end
            if GP.sign > 0, y = y + lgampdf(GP.sign,GP.opts.Psign.alpha,GP.opts.Psign.beta);  end
            if GP.l > 0,    y = y + lgampdf(GP.l,GP.opts.Pl.alpha,GP.opts.Pl.beta);     end   
        end
        function y = sample(GP)
            y = imvnrnd(GP.mu, GP.S_inv);
        end
        function [GP,X] = MCMC(GP,D, nsamples)            
            %% make basic mcmc RW sampler.
            if nargin < 3, nsamples = 10; end
            x0 = log(GP.GP2vec());
            I = isfinite(x0);
              
            J = [GP.opts.Tl, GP.opts.Tsigf, GP.opts.Tsign] >0;
            J = J(I);
            proprnd = @(x)( x + (randn(size(x)) .* J )/10 );
            
            
            
            
            function y = mylogp(lx)
                GP2 = GP; 
                dx = x0;
                dx(I) = lx;
%                 disp('mlogpy');
%                 exp(dx)
                
                GP2 = GP2.vec2GP(exp(dx),false); 
%                 disp('sz');
%                 size(D)
%                 x0
%                 dx
%                 GP.logp(D)
%                 exp(dx)
%                 GP2
%                 keyboard
%                 GP2.logp(D)
                
                y = GP2.logp(D) + sum(lx(isfinite(lx))); 
            end
%             x0(I)
%             x0(I)
            X = mhsample_fast(x0(I),nsamples,'logpdf',@mylogp,'proprnd',proprnd,'symmetric',1);
            x0(end,I) = X(end,:);
            
            GP = GP.vec2GP(exp(x0));
            GP.S_inv = inv(GP.S);
            GP.stats.l(end+1,1) = GP.l;
            GP.stats.sigf(end+1,1) = GP.sigf;
            GP.stats.sign(end+1,1) = GP.sign;
            
        end
    end
end