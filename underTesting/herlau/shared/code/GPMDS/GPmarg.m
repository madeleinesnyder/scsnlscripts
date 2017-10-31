function [MU,S,muProj] = GPmarg(GP,s,I)
% produce marginal of the GP at interval I (horz). 
% the observations s are assumed to be M \times T (standard format).
if false,
    clc
    N = 10;
    M = rand(N); 
    M = M * M'; 
    [L,~] = cholcov(M);
    Mr = M(1:end-1,1:end-1);
    % B1 = L(1:end-1`);
    c = L(end,1:end-1);

    % c = L(1:end-1,1);
    l = L(1:end-1,1:end-1);

    l' * l + c' * c - Mr
    l' * l - Mr
end

%%
if nargin < 1,     
    dT = 1/4; 
    Tscale = 40; 
    M = 1; 
    T = 400; 
    sigf = 1; 
    sign = .1; 
    
    GP = mk_GP(M,T,dT, Tscale, sign, sigf);  % not sure this has to depend on M...
    close all;
    s = GP_rnd(GP);
    M = 3; 
    s = zeros(M,T); 
    for m=1:M, 
        s(m,:) = GP_rnd(GP);
    end
    plot(s'); 
    I = 5;       
    tic()
    for t=1:400,
        [MU,S] = GPmarg(GP,s,I);
    end
    toc();
    if nargout == 0, 
        close all;
        plot(s'); hold on;
        x = mvnrnd(MU,S)
        plot(I,x,'k.-')
    end
end
% compute marginals at $I$.

nStar = setdiff(1:size(s,2),I);

kstar = K(GP, I, nStar); %s(~I,:));    
Kss = K(GP, I);
    
KK = K(GP,nStar); %s(~I,:));
L = cholcov(KK );
%%   
v = L' \ kstar;
S = Kss -v' * v;
      
al2 = L \ (L' \ s(:,nStar)');
MU = GP.mu(I,:)' + (kstar' * al2)';
 
muProj = (kstar' / (L' * L))';


end