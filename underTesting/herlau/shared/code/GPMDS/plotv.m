
function plotv(v,xs)
yl = ylim();
if nargin < 2, xs = 1:length(v); end
% rectangle('Position',[x,y,w,h])
v = mean(v,1);
v = v - min(v(:));

while true, 
    
    j1 = find(v, 1, 'first');
    if isempty(j1), break; end
    v(1:j1) = 1; 
    j2 = find(v==0, 1,'first');
    if isempty(j2), j2 = size(v,2); end
    dH = yl(2)-yl(1); 
    dd =dH/100; 
    dH = dH-2*dd;
    rectangle('Position',[xs(j1)-1,yl(1)+dd, xs(j2)-xs(j1), dH],'FaceColor',.8+[0 0 .1],'EdgeColor',[1 1 1]);  
    v(1:j2) = 0; 
end



end