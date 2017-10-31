function VBslide()
clc; close all;
x = linspace(-5,30,1000)-10;
mu = [0, 2, 4]*5 - 10;
sigma = [1, 2, 2];
K = length(mu);
y = x*0; 
for j=1:K,
    y = y + normpdf(x, mu(j), sigma(j))/K;
end
f = @(z)-sum(normpdf(x,z(1),z(2)) .* log( y ./ normpdf(x,z(1),z(2)) ));
%
M0 = [0, 10, 20]-10;
addpath('../../matlabfrag/');

ff = figure();
plot(x,y,'-','LineWidth',2);box off;
legend({'True distribution $p(x)$'},'Interpreter','latex')
mlf2pdf(ff,sprintf('../../IntroPresentation/figA'));
for g=1:3,
    hold off;
    ff = figure();
    plot(x,y,'-','LineWidth',2);
    [z,val] = fminsearch(f,[M0(g) ; 2]);
    hold on;
    y2 = normpdf(x,z(1), z(2));
    hold all;
    plot(x,y2 * interp1(x,y,M0(g))/max(y2),'-','LineWidth',2);
    box off;
    title(sprintf('KL(p,q) $\\propto$ %2.3f',val/100));
    if g == 2, 
        legend({'True distribution $p(x)$','VB: $q(x)$'},'Interpreter','latex');
    end
    mlf2pdf(ff,sprintf('../../IntroPresentation/figA%i',g));    
    
end
%%
close all;
ff = figure();
xx = linspace(min(x),max(x),40);
S = 300; 
v = xx*0;
for j=1:K,
    [dv,~] = histc(normrnd(mu(j),sigma(j),[1,S]), xx);
    v = v + dv;
end
lp = [0 0];
lp(2) = bar(xx + (xx(2)-xx(1))/2,v * trapz(x,y) / trapz(xx,v),'FaceColor','r' ); hold all;
lp(1) = plot(x,y,'-','LineWidth',2);
xlim([-5,30]-10);
%%
box off;
legend(lp,{'True distribution $p(x)$','MCMC'},'Interpreter','latex');

mlf2pdf(ff,sprintf('../../IntroPresentation/figAA',g));   

end

