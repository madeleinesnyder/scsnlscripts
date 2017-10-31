function presentation_plots()
addpath('../matlabfrag'); 
close all; 
f = figure('Units','Centimeters','Position', [4, 4, 12, 10]);
rng(7);
n = 8;
tt = trnd(1, n, 2);
cc = get(gca,'ColorOrder') ; 
hold on;
pl = []; 
for j=1:size(tt,2)
    dpl = plot(j, tt(:,j)', 'ko', 'MarkerFaceColor',cc(j,:),'MarkerSize',5); 
    pl(j) = dpl(1);
    
    hold all;
end
ylim([-1,1]*10);
grid on;
box off;
for s=1:2,
    lgs{s} =  sprintf('Estimate of $\\beta_{k=%i}$ for each subject', s);    
end
ylabel('Value of $\beta_k$ as estimated for each subject','Interpreter','none'); 
legend(pl, lgs,'Interpreter','latex'); 

mlf2pdf(f, 'plotA'); 
close all;




%%
end