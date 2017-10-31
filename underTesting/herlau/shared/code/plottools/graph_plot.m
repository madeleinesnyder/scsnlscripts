function graph_plot(p,C_m, B_m,opts)
%% plot the matrix C_m and B_m using graph plotting utilities. 
if nargin < 3
    J= 2; 
    M = 5; 
    rng(3); 
    for m=1:M
        p.roi_names{1,m} = sprintf('ROI%i', m);
    end
    C_m = rand(M,M,J)-0.5; 
    B_m = rand(M,M,J) < .5; 
    clc;  close all;
    
    for j=1:size(C_m,3)
        subplot(1,J,j);

        dplot(p,C_m(:,:,j), B_m(:,:,j)); 
    end    
    return;
end
if nargin < 4, opts = struct(); end
%do.print = '
%%k,k,k,,k,

%J = size(C_m,3);
%for j=1:size(C_m,3)
%    subplot(1,J,j);

    dplot(p,C_m, B_m); 
%end

end
function dplot(p, C_m, B_m)
clc
M = size(C_m,1);               
G = digraph();
cc =1;
%if any(J),
    %
    G = addnode(G,p.roi_names);
%end

nn = []; 
for i=1:M
    for j=1:M
        if B_m(i,j) && i ~= j
            G = addedge(G,i, j, C_m(i,j) ); 
            
            nn(cc,1) = i; 
            nn(cc,2) = j;         
            nn(cc,3) = C_m(i,j);        
            cc = cc+1;
        end
    end
end 
 J = find(sum(B_m,1)' + sum(B_m,2) == 0);

 
h = plot(G, 'Layout', 'circle');
%dds = 0.1; 
%h.XData = h.XData*dds;
%h.YData = h.YData*dds;

for j=1:length(p.roi_names), dln{j} = ' '; end
labelnode(h,1:M,dln);

cols = colormap('jet');
for i=1:size(nn,1)
    N = size(cols,1);    
    minmax = [-1,1];
    k = interp1(minmax,[1,N],nn(i,3));
    k = round(k);     
    highlight(h,nn(i,1), nn(i,2),'EdgeColor',cols(k,:),'LineWidth',abs(nn(i,3))*10 );    
end
axis equal;
box off; 
set(gca,'visible','off');     

% h = plot(g,'NodeLabel',[])
hold off;
N = length(h.NodeLabel);
if N == 5, 
    angles = 360*(( [1:N]-1)/N)  + ( [1:N] > 2 & [1:N] < 5)*180;
    angles = (1:N)*0; 
    dy = []; dx = [];
    for i=1:N,
       x = h.XData(i);
       y = h.YData(i);
       dd1 = 0.4;    
       dd2 = 0.42; 
       dx(1,i) = 0;% + (dd1 + (x<0)*dd2)*x;
       %dy(1,i) = (dd1 + (x<0)*dd2)*y; 
       dy(i) = 0;
       if i == 3 || i == 4
            dy(i) = -h.YData(i) + h.YData(i + (-1)^(i == 3)); %dd1 * (-1)^(i==4); 
            dx(i) = dx(i) - 0.6
       end
       if i ~= 3 && i ~= 4, 
           dx(i) = dx(i) + 0.2;
       end
       
       
    end
    
    
else
    angles = 360*(( [1:N]-1)/N)  + ( [1:N] > 2 & [1:N] < 5)*180;    
    for i=1:N,
       x = h.XData(i);
       y = h.YData(i);
       dd1 = 0.12;    
       dd2 = 0.42; 
       dx(1,i) = 0 + (dd1 + (x<0)*dd2)*x;
       dy(1,i) = (dd1 + (x<0)*dd2)*y; 
    end
end


for i=1:length(h.NodeLabel)
   x = h.XData(i);
   y = h.YData(i);
   dd1 = 0.12;    
   dd2 = 0.42; 
   
   text(h.XData(i)+dx(i),y + dy(i),p.roi_names{i},'fontsize',16, 'rotation', angles(i) );
   %text(h.XData(i)+0 + (dd1 + (x<0)*dd2)*x,y + (dd1 + (x<0)*dd2)*y,p.roi_names{i},'fontsize',16, 'rotation', angles(i) );   
end
h.MarkerSize = 12; 
%% 
%set(h,'FontSize',10)
%ar = findall(gca, 'Type','text');
%set(findall(gca, 'Type','text'), 'FontSize', 22)
%ar
%%
%matlabfrag('graph1')
%system('pdflatex graph1')
%mlf2pdf(gcf,'graph1A');


%system
%%
end