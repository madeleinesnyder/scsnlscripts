close all
addpath('matlabfrag');
%addpath('../exportfig');
clc
fpl=figure('units','centimeters','Position',[10 10 10 5]); % [llx lly width height]
xs = linspace(0,2*pi);
plot(xs, sin(xs), 'k-');
hold on; 
M = 200;
zs = randsample(1:length(xs), M, true);

plot(xs(zs), sin(xs(zs)) + randn(1,length(zs))*0.5,'r.');
xlabel('Time $t$/seconds')
ylabel('Signal');
box off
legend('$\sin(t)$', 'Noise');
title('PDF version');
mlf2pdf(fpl,'../../Latex/figures/mout/mplotpdf','SIUnits');

%%
export_fig('../../Latex/figures/mout/mplotpng.png');
title('As eps');
export_fig('../../Latex/figures/mout/mploteps.eps');
3
%% 2ldsfa
