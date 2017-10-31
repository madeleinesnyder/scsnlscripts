function betatest()
clc
close all;
dts = linspace(0.4,2.1); 
dt = [.49,.72,2]; 
B0 = [2.5,3]; 
KC = @(dT,C_cut) ( 1-(1 + dT .* C_cut).^2 );

f = @(beta0,dT,tau)1./(beta0.^2 .* (1-exp(-2*dT./tau))); 
f = @(beta0,dT,tau)1./(beta0.^2 .* KC(dT,tau) ); 
%%
tau_c = -0.2;
dT = 0.72; 
%KC = (1-exp(-2*dT./tau_c));
%clc
%1/(3^2 * KC)^2
%1/(2 * KC)
%1.4^2

%%
c = {'r','b'};
for i=1:length(B0)
    plot(dts,f(B0(i),dts,tau_c),[c{i}, '-']); hold on;        
    plot(dt,f(B0(i),dt,tau_c),[c{i}, '.']); hold on;    
end
y0 = [f(B0(2),dt(1:2),tau_c), f(B0(1),dt(3),tau_c)];
plot(dt,y0,'ko');

function [E,y] = Efun(x)
    beta2 = x(1); tau_c2 = x(2);
    y = f(beta2,dt,tau_c2);
    E= sum((y - y0).^2); 
end
x0 = [B0(2) ; tau_c];
x1 = fminsearch(@Efun,x0);

yy = f(x1(1), dts, x1(2));
plot(dts, yy, 'g-');
T = 200000; 
%%
for i=1:T, 
    x = x0 + 4*randn(2,1);
    x(2) = tau_c;
    xs(:,i) = x;
    EE(1,i) = Efun(x);
end
[E,j] = min(EE); x2 = xs(:,j)
yy = f(x2(1), dts, x2(2));
plot(dts, yy, 'm-');
disp([B0(2), x1(1) x2(1); tau_c, x1(2), x2(2)]);
%%
234;
end