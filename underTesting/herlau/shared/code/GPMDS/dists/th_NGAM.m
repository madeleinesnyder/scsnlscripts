% this is a prior of 
% x ~ N(mu0,sigma^2)     
% and either 
% sigma^2 ~ Gam(alpha,beta)
% or 
% sigma ~ Gam(alpha,beta)

classdef th_NGAM < th_MVN_constrained
    
    properties 
        S_inv
        S
        L
        type
        %alpha 
        %beta
        th_sigma_GAM % this takes care of the sigma or sigma^2 distributions. 
        x
        s0
        tau
        %% new try
        stats
        %% do higher-order things on alpha,beta. 
%         Palpha
%         Pbeta        
    end
    methods
        
        
        function p = th_NGAM(mu0,alpha,beta,opts) % how does the scaling work?
            if nargin < 4, 
                opts = struct();
            end
            
            od.IV = zeros(0,2);
            od.name = 'Normal-Gamma'; 
            
            od.sigma2_on = false;
            od.Tsigma = 1; 
            
            od.Palpha.alpha = 1; % these are nesecary. 
            od.Palpha.beta = 1;                    
            od.Pbeta.alpha = 1; 
            od.Pbeta.beta = 1; 
            %od.Pbeta.sample = false;
            %od.Pbeta.sample = false;
            
            opts = ssfr(od,opts);
            
            
            sigma = ones(1,length(mu0)) * (alpha / beta).^( 1/( od.sigma2_on+1) ) ;
            p@th_MVN_constrained(mu0,sigma,opts);
            
            p.th_sigma_GAM = th_GAM(alpha,beta,opts);
            
            %245
            %% continue with yer stuff here.            
            
            if nargin < 1, 
                addpath('../');
                close all;
                %%
                figure(3);
                Meq = 2;
                Min = 2; 
                dx = 10; 

                M = Min + Meq;
                Sig0 = eye(Min)*2 + (1-eye(Min))*0.8;
                for k=1:Meq,
                    Sig0(Min+k,Min+k) = 1; 
                end
                
                IV = [0,10]; % ; -10, 0 ];
                if Min > 1, 
                    IV(2,:) = [-10,0];
                end
                for j=M-Meq+1:M,
                    IV(j,:) = 2; 
                end
                IV = IV + dx;
                mu = zeros(M,1) + dx; 
                N = 10000;
                opts.IV = IV; 
                p2 = th_NGAM(mu, 1, 3,opts);
                tic()
                [Qx,qx,Rx,rx] = p2.intervals2Qq(p2.IV);
                for i=1:N,
                    while true,  
                        x = mvnrnd(mu(1:Min,:),Sig0(1:Min,1:Min));                        
                        x(M-Meq+1:M) = IV(M-Meq+1:end,1);
                        
                        if all(Qx' * x' - qx < 0), break; end
                    end
                    XA(i,:) = x; 
                end
                toc()
                tic()
                %
                T = 20000; 
                S_inv0 = inv(Sig0);
                XB = zeros(T,M); 
                XB(1,:) = x; 
                for t=2:T,
                    XB(t,:) = p2.post_sample(XB(t-1,:)',S_inv0, Sig0 \ mu);                
                end
                toc()
                hold on;
                if Min == 2, 
                plot(XA(:,1),XA(:,2), 'k.');                                
                plot(XB(:,1),XB(:,2), 'r.');
                end
                
                [ mean(XA,1) ; mean(XB,1) ] - dx
                %% 
                % Do second experiment. Does 4 samplers agree?
                return;
                figure(1);   
                K = 20; 
                N = 100;
                T = 200; 
                
                for k=1:K,
                    IVSUP(k,:) = [0,1] - (rand()<.5);                    
                end
                
                opts.IV = IVSUP;                
                p2 = th_NGAM(zeros(K,1), 1, 3,opts);
                
                [~,X] = p2.sample(N);

%                 close all;
                 
                plot(1:K,X,'k.');
                %%      
                figure(2);
                lgs = {'X'}; c= 1;
                sig0 = sqrt(mean(X.^2,1)); %,1,1);
                plot(1:K, log(sig0),'ko'); hold all;
                cmap = get(gca,'ColorOrder');
                plts = zeros(1,4);
                for sig2_on = 1:2,
                    for IV_on = 1:2,  
                        
                        aa = 0.3; 
                        opts = struct();
                        if IV_on == 2,  
                            opts.IV = IVSUP;
                        end
                        
                        opts.sigma2_on = sig2_on == 2;
                        p2 = th_NGAM(zeros(K,1), .001, 1,opts);
                        
                        for t=1:T,
                            ss(t,:) = p2.sigma;
                            p2 = p2.MCMC(X);
                        end 
                        ss = ss(end/2:end,:);
                        lgs{c} = sprintf('sig2=%i, IV=%i',opts.sigma2_on, size(p2.IV,1));
                                                
                        px = plot( (1:K) + (c-1)/6, log(ss),'.','Color',cmap(c,:)); hold all;
                        px = plot( (1:K) + (c-1)/6, log(mean(ss,1)),'*','Color',cmap(c,:)); hold all;
                        
                        plts(c) = px(1); %plot( (1:K) + (c-1)/6, log(ss),'.','Color',cmap(c,:))
                        c = c+1;
                    end
                end
                legend(plts,lgs);                
                return; 
            end
             
            p.M = size(mu0,1);            
            p = p.set( ones(1,length(mu0))*p.th_sigma_GAM.alpha/p.th_sigma_GAM.beta);
                                    
            %p.alpha = alpha;
            %p.beta = beta;
            p.type = 7;
            p.stats.sigma = zeros(1,p.M);
            s1 = '';
            s2 = '';
            if opts.sigma2_on, s1 = '^2'; end
            if size(p.IV,1) > 0, s2 = ', (Box-restricted)'; end
            
            p.name = sprintf('Normal-Gamma, %s, sigma%s ~ Gam(alpha,beta)%s',opts.name,s1,s2);
            disp(p.name)
            p.mu = mu0;
        end
        function p = set(p,sigma)
            p.sigma = sigma;
            p.S_inv = diag(1./p.sigma.^2);
            p.S = diag(p.sigma.^2);
            p.L = cholcov(p.S);
        end        
        function [p,x] = sample(p,N)    
            if nargin < 2, N = 1; end
            if p.opts.sigma2_on, 
                sigma = gamrnd(p.th_sigma.GAM.alpha,1./p.th_sigma.GAM.beta,[1,p.M]) .^ .5;
            else
                sigma = gamrnd(p.th_sigma.GAM.alpha,1./p.th_sigma.GAM.beta,[1,p.M]);                
            end
            
            p = p.set(sigma);  
            x = zeros(1,p.M);
            for i=1:N,
                for k=1:p.M,
                    while true, 
                        z = normrnd(p.mu(k)',p.sigma(k));
                        if size(p.IV,1) == 0 || (p.IV(k,1) <= z && p.IV(k,2) >= z), break; end
                    end                
                    x(i,k) = z; 
                end              
            end                        
        end  
        function [y,yvec] = logp(p,X)
            %%
            yvec = 0; 
            I = true(p.M,1);
            if size(p.IV,1) > 0,
                I = p.IV(:,1) ~= p.IV(:,2); 
            end 
            
            sig = p.sigma(:,I);    
            yvec = zeros(1,p.M);
            yvec(I) = sum(lnormpdf(X(:,I),p.mu(I,:),sig'),1);
            %%
            pp = 1 + (p.opts.sigma2_on);      
            [~,lpGvec,lpGrest] = p.th_sigma_GAM.logp(sig'.^pp);
            
            yvec(I) = yvec(I) + lpGvec'; %lgampdf(sig .^ pp,p.alpha,p.beta);            
            if size(p.IV,1) > 0 && false, % this is not good
                assert(false);
                sigma = p.sigma2.^.5;  
                lcdf = normcdf(p.IV(:,1), p.mu, sigma');
                ucdf = normcdf(p.IV(:,2), p.mu, sigma');
                
                yvec = yvec - log(ucdf - lcdf)';
                 
            end
            y = sum(yvec) + lpGrest;            
            
        end
        
        %% alpha/beta needs to be sampled. 
        
        
        function X2 = post_sample_DEFUNCT(p,X,S_inv_all,q_all)
            %%
%             load('dbug_NGAM')
            
            % make posterior sampling of X. 
            if size(p.IV,1) > 0, 
                [Qx,qx] = p.intervals2Qq(p.opts.IV);
                %% Mikkel paper contains an inaccuracy:
%                 Sig = inv(S_inv_all);
%                 L = cholcov(Sig)'; 
                % z  = (L')\(x-mu);
                %%
                J = eye(size(S_inv_all,1));
                J = J(:,end:-1:1);                                        
                JAJ = J * S_inv_all * J;                    
                %%
                L2 = cholcov(S_inv_all);
%                 L2' * L2 - S_inv_all
                %%
                LL = cholcov(JAJ)';
                LL = (J * (LL\J))';
                    
                x = X(:);
                mu = S_inv_all \ q_all; 
                %z = (LL')\(x-mu);
                z  = L2 * (x-mu);
                
                if false,
                    clc
                    xmu = x-mu;
                    z1 = (L') \ xmu; 
                                        
                    J = eye(size(S_inv_all,1));
                    J = J(:,end:-1:1);                                        
                    JAJ = J * S_inv_all * J;                    
                    LL = cholcov(JAJ)';
                    LL = (J * (LL\J))'; 
                    LL * LL' - Sig
                    
                    L * L' - Sig
                    
                    %%
                end 
                
%                 Qz = LL * Q;
                Qz = (L2') \ Qx;
                qz = qx - Qx' * mu;
                for i=1:p.M,
                    d = Qz(i,:)';
                    n = qz - Qz(1:end ~= i,:)' * (z(1:end ~= i,:) );
                    nd = (n ./ d);
                    l = max( [-inf ; nd(isfinite(nd) & d < 0) ] );
                    u = min( [ inf ; nd(isfinite(nd) & d > 0) ] );
                    
                    lu(i,:) = [l,u];
                    luinv = normcdf([l,u],0,1);
                    if luinv(1) ~= 1 && luinv(2) ~= 0, 
                        z(i) = norminv( rand() * (luinv(2)-luinv(1) ) + luinv(1) );                
                    end
                    if luinv(1) > 1-10^-10
                        z(i) = l;
                    end
                    if luinv(2) <10^-10
                        z(i) = u; 
                    end                    
                    if ~isfinite(z(i))
                        234
                    end 
                end
                X2 = L2\z + mu;
                %%
            else
                X2 = imvnrnd(S_inv_all \ q_all, S_inv_all);                
            end
            %%
            % load('dbug_NGAM')
        end
        
        function p = MCMC(p,X,SAMPLES)
            if nargin < 3, SAMPLES = 4; end
            %% do either exact sampling or random walk stuff. 
            sigma_initial = p.sigma; 
            %p.name
            p.th_sigma_GAM = p.th_sigma_GAM.MCMC(p.sigma');
             
            if p.opts.sigma2_on, 
                [N] = size(X,1);
                SS = sum(bsxfun(@minus, X, p.mu').^2,1);
                sig2 = p.sigma.^2;
                for i=1:p.M,
%                     if i == 20,
%                         234
%                     end
                    if size(p.IV,1) > 0 && p.IV(i,1) == p.IV(i,2), continue; end
                    gi_p = p.th_sigma_GAM.alpha - N/2;
                    gi_a = 2*p.th_sigma_GAM.beta;
                    gi_b = SS(i);
                    if SS(i) == 0,
                        sig2(i) = gamrnd(p.th_sigma_GAM.alpha,1./p.th_sigma_GAM.beta);                    
                    else
                        sig2(i) = gigrnd(gi_p, gi_a, gi_b, 1);                
                    end
                    if sig2(i) == 0 || isnan(sig2(i)), 
                        assert(false);
                    end
                end
                p = p.set(sig2.^.5);            
            else 
                %%
                for j=1:SAMPLES, 
%                     [j,size(X)]
                    if p.opts.Tsigma > 0, 
                        if size(p.IV,1)>0,                            
                            I = p.IV(:,1) ~= p.IV(:,2);
                        else
                            I = true(p.M,1);
                        end
                        lsig_0 = log(p.sigma);
                        lsig_1 = lsig_0 + randn(size(lsig_0))/2;
                        [~,yv_0] = p.logp(X);
                        p1 = p;
                        p1 = p1.set(exp(lsig_1)); 
                        [~,yv_1] = p1.logp(X);
                        acc = min(1,exp(yv_1 - yv_0 + lsig_1 - lsig_0));
                        Iacc = rand(size(acc)) - acc < 0;
                        lsig_0(Iacc & I') = lsig_1(Iacc & I');
                        p = p.set(exp(lsig_0)); 
                    end  
                    %%
                    % make image of the marginals. 
                    if false, 
                        %% 
                        close all;
                        T= 500;
%                         ss = linspace(0.001,.1, 1000);
                        sx = zeros(p.M,T);
                        
                        y = []; 
                        y2 = [];
                        for k=1:length(p.sigma),
                            sx(k,:) = linspace(std(X(:,k))/4, std(X(:,k))*4, T);
                            for i=1:T,
                                s2 = p.sigma;
                                p2 = p; 
                                s2(k) = sx(k,i);
                                p2 = p2.set(s2);
                                y(i,k) = p2.logp(X);         
                                yB(i,k) = sum(lnormpdf(X',0,sx(k,i)) ) - log(sx(k,i));
                            end                
                        end
%                         close all;
    %                     y = exp(y);
                        for k=1:p.M,
                            H = ceil(sqrt(p.M)); W = ceil(p.M/H);
                            subplot(H,W,k);
                        
                            y2 = y(:,k); y2= exp(y2 - max(y2)); 
                            yB2 = yB(:,k); yB2 = exp(yB2-max(yB2));
                            plot(sx(k,:), y2); hold all;
                            plot(sx(k,:), yB2,'.'); hold all;

                            plot(p.sigma(k), 0, 'ko');
                            plot(sqrt(mean(X(:,k).^2)), 0, 'r*');
%                             plot(0,hh,'.');

                        end
                    end
                    %%
                end
                
            end
            %% reject update if less than threshold value. 
            Sig_min = 10^-6; 
            sigma = p.sigma;            
            sigma(sigma < Sig_min) = sigma_initial(sigma < Sig_min);
            sigma = max(Sig_min,sigma); % just in case the initailization is bad. 
            p = p.set(sigma);
            
            if ~isfield(p.stats,'sigma'),
                p.stats.sigma = zeros(size(p.sigma2));
                       
            end
            if ~isfield(p.stats,'sigma_T'), 
                 p.stats.sigma_T = 0;        
            end
             
            p.stats.sigma = (p.stats.sigma * (p.stats.sigma_T) +  p.sigma) / (p.stats.sigma_T+1); 
            p.stats.sigma_T = p.stats.sigma_T+1;
            
        end        
    end
end