
function KendallMatrix()
% table1();
table2();
end
function table2() % do J=2 results here. 
s(1).dir = '../hcp/hcp_J2_processed_A1';
for i=1:length(s),
    s(i).x = load(s(i).dir);
end
% what should the table look like?
MODE = 1; % mode (1==MAX, 2=ij. 
J = size(s(1).x.D{1}.C_m{1},3);

for j=1:length(s),
    %break
    sz_sims = size(s(j).x.D);    
    s(j).Km_cv = zeros( sz_sims );
    s(j).Ksd_cv = zeros( sz_sims );
    s(j).Kmax = zeros( sz_sims );
    
    for m=1:numel(s(j).x.D),
        [mx,my] = ind2sub(size(s(j).x.D),m);
        D = s(j).x.D{m};
        Cms = D.C_m;
        
        subs = size(D(1).C_m,1);
        Kfolds = 3; 
        II = crossvalind('Kfold',subs,Kfolds);
        
        R = zeros([size(Cms,2),size(Cms,3) 2]);
        globalK = zeros(size(Cms,2),size(Cms,3));
        yfolds = zeros(1,Kfolds);
        
        %% instead of doing K-fold CV we will do bootstrap.
        mall = zeros(size(D.C_m,2), size(D.C_m,3), J);
        for ii=1:size(D.C_m,2),
            for jj=1:size(D.C_m,3),
                [~,~,mall(ii,jj,:)] = matrix_consistency(D.C_m(:,ii,jj));
            end
        end
        
        %%
        [~,k] = max(vec(mall(:,:,1)));
        [ii,jj] = ind2sub(size(mall(:,:,1)), k);
        %%
%         [ymax,yboot,mall] = Kestimate(Cms);
        s(j).sim(mx,my).Kmax= mall(ii,jj,:);
        
%        s(j).sim(mx,my).Kboot = yboot;
        s(j).sim(mx,my).mall = mall;
        
        s(j).sim(mx,my).ij = [ii,jj];

        ii=2;jj=2;
        s(j).sim(mx,my).Kijg = mall(min(size(mall,1),ii),min(size(mall,2), jj),:);         

    end
end
%%
NAMES = {}; for i=1:length(s), NAMES{i} = s(i).x.name; end
dSIMS = {'DCM','MDS','MDS-MCMC'}; 
SIMS = {}; 
M = {};
B = {};
fld = 'Kmax';
for j=1:J,
    SIMS = {SIMS{:}, dSIMS{:}};
end
M = zeros(1,length(SIMS)); 
B = zeros(size(M));
for i=1:length(s),
    dM = {}; 
    v = permute(cell2mat({(s(i).sim(1,:).(fld)) }) , [2, 3, 1]);
    
    [~,dx] = max(v,[],1);
    dM = v;
    dB = v*0; 
    for q=1:J,
        dB(dx(q),q) = true; 
    end
    M(i,:) = v(:)';
    B(i,:) = dB(:);
    %end
    %cell2mat(dM)    
end
[H,W] = size(M);
M(end+1,:) = 0;
B(end+1,:) = 0;
%
NAMES{end+1} = 'LAST';
save('kendall_tbl2', 'NAMES','SIMS','M','B','H','W');
matlab_jinjafy('kendallJ2_partial.tex');

end
function table1(), % this table represent the J=1 results.
%% load datasets of various kinds.
s(1).dir = '../oddball/oddball_processed_A1';
s(2).dir = '../Yale/yale_processed_A1';
s(3).dir = '../StanfordStopGo/stopgo_processed_A1';
s(4).dir = '../hcp/hcp_processed_A1';

for i=1:length(s),
    s(i).x = load(s(i).dir);
end
%%
Kfolds = 3;
for j=1:length(s),
    sz_sims = size(s(j).x.D);
    
    s(j).Km_cv = zeros( sz_sims );
    s(j).Ksd_cv = zeros( sz_sims );
    s(j).Kmax = zeros( sz_sims );
    
    for m=1:numel(s(j).x.D),
        [mx,my] = ind2sub(size(s(j).x.D),m);
        D = s(j).x.D{m};
        Cms = D.C_m;
        
        subs = size(D(1).C_m,1);
        II = crossvalind('Kfold',subs,Kfolds);
        
        R = zeros([size(Cms,2),size(Cms,3) 2]);
        globalK = zeros(size(Cms,2),size(Cms,3));
        yfolds = zeros(1,Kfolds);
        
        %% instead of doing K-fold CV we will do bootstrap.
            
        [ymax,yboot,mall] = Kestimate(Cms);
        s(j).sim(mx,my).Kmax= ymax;
        s(j).sim(mx,my).Kboot = yboot;
        s(j).sim(mx,my).mall = mall;
        s(j).sim(mx,my).Kijg = mall(1,min(size(mall,2), 2));         
        
    end
end
%%    
% make results table. How should it be arranged? As a per-K case thing?
fields = {'Kmax', 'Kijg','Kboot'};
for fi=1:1,%length(fields),
    ss = '';
%     field = fields{fi};
    Ms = [2,3,5];
    for i=1:length(s),
        for kk=2:3,
            
            %sA = [sprintf('& %1.2g ', [s(i).sim(kk,:).(field)]), '\\\\'];
            %%
%             s(i).
%             boldv(v)
        sF = {};
        for q=1:length(fields),
            sF{q} = boldv([s(i).sim(kk,:).(fields{q})]);
        end       
            sF = [sprintf('%s ', sF{:}), '\\\\'];            
            ss = sprintf('%s %s ($M=%i$) %s\n',ss, s(i).x.name, Ms(kk), sF);            
        end
    end
    f = fopen(sprintf('../Latex/tbl_Kendall'),'w');
    fprintf(f,ss);
    fclose(f);
end
disp(ss)

%%
% Make another plot. 
close all; 
f = figure('Units','Centimeters','Position',[2,2,36,10]);
Ks = [2,3]; cc=1; 
for k=1:length(Ks),
    for i=1:length(s),
        subplot(length(Ks),length(s), cc);
        m = s(i).sim(Ks(k),3).mall; 
        imagesc( m );
        axis image;
        title(sprintf('%s,K=%i',s(i).x.name, Ks(k) ));
        colorbar;
        cc = cc+1;
    end
end
doprint = true;
if doprint,
    addpath('../../matlabfrag/');
    mlf2pdf(f, '../Latex/figures/KendallAll');
end

%%
end
function s = boldv(v)
    [~,j] = max(v);
    s = '';
    for i=1:length(v)
        ds = sprintf('%1.2g ', v(i));
        if i == j,
            ds = sprintf('\\\\textbf{%1.2g} ', v(i));

        end
        s = [s, '& ' ds];
    end
end
% format: 
% Cms = cell : Subjects x cond1 x cond2
% goal is to find max condition using some cv like procedure.
function [ymax,yboot,mall] = Kestimate(Cms)
%% build consistency matrix pairs of all elements in the cell array.

% Cms = Cms

% Cms = Cms(:,1,:);
sz = size(Cms);
if length(sz)<3,sz(3)=1; end
KMATS = cell(sz(2:end));

% if sz(2) == 1, 
    
   
% else
   for i=1:sz(2),
       for j=1:sz(3),
           [~,KMATS{i,j}] = matrix_consistency(Cms(:,i,j));
       end
   end
   if prod(sz(2:end)) == 1,
        ymax = utrisum(KMATS{1});
        yboot = ymax;
        mall = ymax;
        return;
   end

   
% end
%% 
% tic()
T = 2000; 
k = 2; 
S = sz(1);
y = zeros(1,T);
ij = zeros(2,T);
mm2 = zeros(3,T);
for t=1:T,
    I = randperm(S);
    Ik = I(end-k+1:end);
    I = I(1:end-k); 
    I = sort(I); Ik = sort(Ik);
    %     S - length(I)
    m0 = cellfun(@(m)utrisum(m,I), KMATS);
    [v,i,j] = max2(m0);
    
    m2 = cellfun(@(m)utrisum(m,Ik), KMATS);
%      i = 1; j = 2;  
    y(1,t) = m2(i,j);
%     disp([v,m0(i,j), y(1,t), i, j])
%     ij(:,t) = [i,j]; 
    
    
%     mm2(:,t) = m2;
end
mall = cellfun(@(m)utrisum(m), KMATS);
[ymax,i,j] = max2(mall);
% [ymax,mean(y)]
yboot = mean(y);
% mean(ij,2)
% toc();
%%

end
function [v,i,j] = max2(m)
    [v,dx] = max(m(:));
    [i,j] = ind2sub(size(m),dx);

end

function y = utrisum(M,I);
if nargin > 1, 
    M = M(I,I);
end
V = triu(ones(size(M)),1);
y = mean(M(V(:)==1));
end