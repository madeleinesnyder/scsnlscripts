function D2 = LDS_plotter(D,opts),
if nargin < 1, 
    yale_DCM_MDS_gridA_plot();
    return;
end
addpath('../../matlabfrag/');
do.Dout = [opts.name, '_LDSprocessed'];
do.doprint = false;
opts = ssfr(do,opts);

% fck the options struct. 
cc = 1; 
MS = {'DCM','MDS','MDS-MC'};
Ms = size(D,1); methods = size(D,2);
for M=1:size(D,1),
    for met = 1:size(D,2),
    figure(1);
    subplot(Ms,3,cc); 
    
if met <=2, 
    v = D{M,met}.C_m;
    Xsd = matSD(v);
    
    [~,~,Kall(M,met,:)] = matrix_consistency(v);
    D{M,met}.K = Kall(M,met); 
    D{M,met}.sel = 1;
    D{M,met}.C_best = D{M,met}.C_m;
else
    cls = {};
    for i=1:size(D{M,met}.C_m, 2),
        for j = 1:size(D{M,met}.C_m, 3),
            Clist = squeeze(D{M,met}.C_m(:,i,j)); 
            
            
            cls{i,j} = Clist;
            [~,~,KK(i,j,:)] = matrix_consistency(Clist);
            MM{i,j} = flatmat(Clist, 0.1, true);
            
        end 
    end
    dK = KK(:,:,1);
    [~,jj] = max(dK(:));
    
    
    D{M,met}.C_best = cls{jj};
    D{M,met}.K = KK;
    for i=1:length(D{M,met}.C_best),
        D{M,met}.C_best{i} = rd(D{M,met}.C_best{i});
    end
        
    end

    end
end
%%
close all;
f1 = figure(1);cc = 1; 
J = size(D{1}.C_best{1},3);

for M = 1:Ms,
    for j=1:J,
    for met =1:methods,
        %%
        Cs = D{M,met}.C_best;
        Cs = reshape(Cs, [1, 1, length(Cs)]);
        [K,~,KJ(met,:)] = matrix_consistency(Cs);
        
        for h=1:length(Cs),
            Cs{h} = Cs{h}(:,:,j);
        end
        
        ds = sprintf('%.2g',KJ(met,j));
        ttl = sprintf('%s (j=%i, %s)', MS{met}, j, ds(1:end));
        
        
%         D{M,met}.K = K;
        figure(1);
        subplot(Ms*J, methods, cc);
        
        
        
        Cs_flat =cellfun(@(m)m(:,:,:),rd(Cs),'UniformOutput',false);
        imagesc(mean(cell2mat(Cs_flat),3)); axis image;
        axis image;
        title(ttl);
        
        f2 = figure(2);
        subplot(Ms*J, methods, cc);
        Xsd = matSD(Cs);
        imagesc(Xsd);
%         title(sprintf('%s, %2.2g', MS{met}, K));
        
        axis image;
        n = size(Cs{1},1);
        %%
        v = log(0.05 / (n*(n-1) ));
        imagesc(flatmat({Xsd, (Xsd > -v) * 8}, 1));
        set(gca,'clim',[0, 8]);
        axis image;
        title(ttl);
        %%
        cc = cc+ 1; 
        %%
        end
    end
end
%% 
if opts.doprint,
    bn = sprintf('figures/%s_f',opts.name);
    pwd
    mlf2pdf(f1,[bn,'1']);
    mlf2pdf(f2,[bn,'2']);
end
figure();
J = size(KK,3); 
for k=1:J, subplot(1,J,k); 
    imagesc(KK(:,:,k)); axis image;
    title('Matrix Consistency'); colorbar;
end
%%
% D2 = cell(size(D));
% for i=1:numel(D),
%     D2{i}.K = D{i}.K;
%     D2{i}.C_best = D{i}.C_best;
    
% end
name = opts.name;
save(opts.Dout,'D', 'name');
disp('Saving processed data to: ');
disp(opts.Dout);
end
function Xsd = matSD(Clist)
XX = cell2mat(reshape(Clist, [1, 1, length(Clist)]));
n = size(XX,1);
K = size(XX,3);
II = repmat(~eye(n), [1, 1, K]);

XX = XX - median(XX(II(:)));
% XX = XX - 1*mean(XX(II(:)));
Xm = mean(XX,3);
Xtilde = sqrt(sum(bsxfun(@minus, XX, Xm).^2, 3)/(K*(K-1)));
nu = K - 1;
alpha = 0.05; 
%%
pX = zeros(size(XX,1));
for i=1:size(XX,1),
    for j=1:size(XX,1),
        [~, pX(i,j)] = ttest(squeeze(XX(i,j,:)), 0, 0.05, 'both');
    end
end
pX = -log(pX);
pX = pX - diag(diag(pX));
Xsd = max(0, pX +0*log(0.01));
end