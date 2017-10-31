function GP = mk_GP(M,T,l,sign, sigf,zero_lock,mu0)
if nargin < 1, test(); return; GPmarg(); return; end
if nargin < 5, sigf = 1; end
if nargin < 6, zero_lock = false; end
if nargin < 7, 
    mu0 = zeros(T,1); 
end
if zero_lock, mu0 = [0 ; mu0 ; 0]; end

% sign = 0.001;
GP = struct();
GP.mu = mu0; %zeros(T + 2*zero_lock,1); % zero mean.
GP.sigf = sigf;
x = (1:T+2*zero_lock)';
%%
GP.l = l; 
GP.sigf = sigf;
GP.sign = sign;
GP.S = K(GP, x);
GP.M = M;
GP.T = T;
GP.zero_lock = zero_lock;
% GP.Tscale = Tscale;
if GP.zero_lock,
%     close all;
    [MU,S] = GPmarg(GP, GP.mu', [2:size(GP.mu,1)-1 ]);
    %% 
    if false,
        for i=1:10, 
            x = [0, mvnrnd(MU,S), 0];
            plot(x); hold all;
        end
        plot([0, MU, 0]);
    end
    %%
    GP.S = S;
    GP.mu = MU';
end

%GP.K = K(GP,(1:T)');
GP.S_inv = inv(GP.S);
if sign == inf,
    GP.S_inv = zeros(size(GP.S_inv));
end

GP.L = cholcov(GP.S);
GP.type = 2; 
GP.name = 'Gaussian Process';
end