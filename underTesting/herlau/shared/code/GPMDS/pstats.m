function pstats(p0,stats)
if nargin < 1, test(); return; end
if ~isfield(stats,'C'), 
    disp('no stats, breaking');
    return; end
%% close all; 
MDS_selfcheck(p0);
f = ifset(1, [1,1,24,20]); 
H = 3; W = 4; 
sf = stats.sf;
spr = struct();
i = 1; 

Cs = reshape(mean(stats.C(:,ceil(end/4*3):end),2),size(p0.C));    
[M,~,J] = size(p0.C);
for j=1:J,    
    subtightplot(H,W,i); 
    s0 = std(p0.s(1).s,1,2)'; 
    s0 = s0/max(s0);
    imagesc([Cs(:,:,j) bsxfun(@times, Cs(:,:,j), s0) ]); 
    colorbar; axis image;
    i = i + 1;    
end

subtightplot(H,W,i); plot(stats.C','-'); i=i+1;
subtightplot(H,W,i); plot(stats.lp(ceil(1):end)); i=i+1;

subtightplot(H,W,i); plot(real(stats.C_eig')); i=i+1; title('eig-real');
subtightplot(H,W,i); plot(imag(stats.C_eig')); i=i+1; title('eig-imag');

subtightplot(H,W,i); imagesc([p0.s.B]); colorbar; i = i+1; title('B');

%%
sU = stats.s(1).U;
sU = sU(sum(sU,2) ~= 0,:);
%%
subtightplot(H,W,i); plot(sU'); i = i+1; title('U');
% subtightplot(H,W,i); imagesc([p0.U]); i = i+1; title('U');
%%
subtightplot(H,W,i); hold off;
cc = get(gca,'ColorOrder');
Sm = min(3,length(p0.s));

for k=1:Sm;
    %%
    hold off;
    v = zeros(size(p0.s(1).Phi));    
    v(1,:) = 1; 
    x = stats.s(1).Phi(v(:)==1,:); 
    m = mean(x,2);
    er = std(x,1,2); 

    %errorbar(m,er,'LineWidth',2); hold all;
    plot(p0.s(k).Phi','k:');  hold on; 
    
%     hold off; plot 
    v = p0.s(k).B * p0.s(k).Phi;
    plot(v','-'); hold all;    
    
end 
    
%%    
i = i + 1; 
%%
subtightplot(H,W,i); hold off;
plot(stats.Pe_mu','-'); hold all;
title('Pe.mu'); 
i = i+1;
%%

%%
% close all;
f2 = ifset(2,[21,1,20,20]);
cc = [cc ; cc ; cc];
for k=1:Sm,
    HH = 6;  
    
    subtightplot(HH,1,1); hold off;
    for j=M:-1:1,
        plot(p0.s(k).s(j,:)','Color',cc(j,:));  hold on;
%         plot(p0.s(k).sf(j,:)',':', 'Color',cc(j,:));  hold on;
    end
    UU = mean(p0.s(k).u(1,:,:),3)';
    plot(UU,'k-');
    
    subtightplot(HH,1,2); hold off;
    for j=M:-1:1,
%         plot(p0.s(k).s(j,:)','Color',cc(j,:));  hold on;
        plot(p0.s(k).sf(j,:)','-', 'Color',cc(j,:));  hold on;
    end
    plot(UU,'k-');
    
    subtightplot(HH,1,3); hold off;
    for j=M:-1:1,
%         plot(p0.s(k).s(j,:)','Color',cc(j,:));  hold on;
        plot(p0.s(k).sf(j,:)' + p0.s(k).s(j,:)','-', 'Color',cc(j,:));  hold on;
    end
    plot(p0.s(k).u(1,:)','k-');
    if p0.Pw_sigma_temporal,
        %%
        hold off;
        plot(mean(p0.Pw(1).stats.sigma( ceil(end/2):end,:),1));         
    end
    
    subtightplot(HH,1,4); hold off;
    for j=M:-1:1, 
        plot(p0.s(k).yhat(j,:)','Color',cc(j,:));  hold on;
        plot(p0.s(k).y(j,:)','.', 'Color',cc(j,:));  hold on;
    end
    plot(UU,'k-');
    plot(p0.s(k).v(2:end,:)','k-'); title('y, yhat');     
    
    subtightplot(HH,1,5); hold off;
    for j=M:-1:1, 
        plot(p0.s(k).w(j,:)','Color',cc(j,:));  hold on;
    end    
    title('noise w'); 
    
    subtightplot(HH,1,6); hold off;
    for j=M:-1:1, 
        plot(p0.s(k).e(j,:)','Color',cc(j,:));  hold on;
    end    
    title('noise e');     
    plot(UU,'k-'); 
end
% i = i -1; 
% figure(f);
%%
% close all;
f3 = ifset(3,[10,1,20,16]);
% set(f3,'Units','Centimeters','Position',[10,1,20,16]);
H=3; W=3; 
c = 1; 
subplot(H,W,c); plot(p0.Pw(1).stats.sigma);  title('Pw.sigma'); c=c+1; 
subplot(H,W,c); plot(stats.Pw_std');  title('Pw.std'); c=c+1; 
%%
subplot(H,W,c); plot(p0.PC(1).stats.sigma);  title('PC.sigma'); c=c+1;
r = ([p0.Psf.stats]); %[sign]
subplot(H,W,c); plot([r.sign]);  title('Psf.sign'); c=c+1;
subplot(H,W,c); plot([r.sigf]);  title('Psf.sigf'); c=c+1;
subplot(H,W,c); plot([r.l]);  title('Psf.l'); c=c+1;

% c = c+1;
r = [p0.Pe.stats];
subplot(H,W,c); plot([r.sign]);  title('Pe.sign'); c=c+1;
subplot(H,W,c); plot([r.sigf]);  title('Pe.sigf'); c=c+1;
subplot(H,W,c); plot([r.l]);  title('Pe.l'); c=c+1;
%%
pause(0.01); drawnow; 

end
function f = ifset(fid,pos)
if ~ishandle(fid),        
    f = figure(fid);
    set(f,'Units','Centimeters','Position',pos);
else 
    f = figure(fid);
end
end