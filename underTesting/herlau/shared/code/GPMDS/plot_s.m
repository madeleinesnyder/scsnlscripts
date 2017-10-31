function plot_s(s,v), 
[M,T] = size(s);
cc = get(gca,'ColorOrder');
cc = [cc ; cc ; cc];
% k = 1; 

plotv(v); hold on;

for j=M:-1:1,
    plot(s(j,:)','Color',cc(j,:));  hold on;
end
    
end
