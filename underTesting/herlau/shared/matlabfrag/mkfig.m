function mkfig()
f = figure('Units','Centimeters','Position',[4, 4, 20, 8]); 
plot(sin(linspace(0,pi)));
xlabel('$\sin(\theta)$');
mlf2pdf(f, 'myfig');
end