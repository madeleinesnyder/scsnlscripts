function texax()
ax = gca;
xticks = get(ax,'XTickLabel');
xtv = get(ax,'XTick');
xticks = fixax(xticks,xtv);
set(ax,'XTickLabel',xticks);

xticks = get(ax,'YTickLabel');
xtv = get(ax,'YTick');
xticks = fixax(xticks,xtv);
set(ax,'YTickLabel',xticks);
% return;
%%
% for i=1:length(xticks),
%     xticks{i} = ['$' xticks{i}  '$'];     
% end
% yticks = get(ax,'YTickLabel');
% for i=1:length(yticks),
%     yticks{i} = ['$' yticks{i}  '$'];     
% end
% set(ax,'YTickLabel',yticks)
end
function xticks = fixax(xticks,xtv)
if any(xticks{1}=='}'),
    for i=1:length(xticks), xticks{i} = sprintf('$%s$',xticks{i}); end
else
     
if mean(abs(cellfun(@(i)str2num(i),xticks)' - xtv)) > 10^(-5),
    %%
    for i=1:length(xtv),
    x = xtv(i);    
    b = floor(log10(abs(x)));
    a = x * 10^(-b);
    s = sprintf('$%g\\cdot 10^{%g}$',a,b);
    if x == 0, s = '$0$'; end
    xticks{i} = s;
    
    end
end
end
%%
end