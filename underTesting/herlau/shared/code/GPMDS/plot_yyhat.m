function plot_yyhat(y,yhat,Pe,v), 
[M,T] = size(y);
cc = get(gca,'ColorOrder');
cc = [cc ; cc ; cc];
plotv(v); hold on;
for j=M:-1:1, 
    if ~isempty(yhat),
    plot(yhat(j,:)' + Pe(j).mu,'Color',cc(j,:),'LineWidth',2);  hold on;
    end
    plot(y(j,:)','-', 'Color',cc(j,:));  hold on;
end 
    
end
