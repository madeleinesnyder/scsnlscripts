function th_GP_test()
clc;
rng(1);
sigf = 5; 
sign = 0.05; 
T = 200; 
mu0 = zeros(T,1);
l = 5;  
gptype = 2; 
GP = th_GP(gptype,mu0,l,sigf,sign,false);
%%
D = GP.sample();  
%%  
GP.logp(D)
x0 = GP.GP2vec();

close all;
subplot(3,1,1);
plot(D,'.-')
%% 
GP.logp(D)
%%
[GP2, X] = GP.MCMC(D,200);
subplot(3,1,2); 
plot(X);
%%
x1 = GP2.GP2vec();
GP2.logp(D)
%%
subplot(3,1,3);
hist(X,100)

[mean(exp(X),1) ; x0]
%%
end