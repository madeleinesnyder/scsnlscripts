function test()
cdir = fileparts(mfilename('fullpath'));
addpath(fullfile(cdir, '../rats'));

% rat();
%%
rng(1)
T = 20; 
M = 3; 
J = 2; 
s = rand(M,T);
C = rand(M,M,J);
v = rand(J,T);
Rb = zeros(T,T);
t = 7;
et = zeros(T,1); et(t) = 1; 
for tt=1:T-1,
    Rb(tt,tt+1) = 1; 
end


Cts1 = zeros(M,1);
for j=1:J,
    Cts1 = Cts1 + C(:,:,j) * s(:,t-1) * v(j,t);
end

Cts2 = Cts1*0;

for j=1:J,   
    V = diag(v(j,:));
    Cts2 = Cts2 + C(:,:,j) *  s * Rb* V * et;
end
% [Cts1, Cts2]

%%
% return;
addpath('../rats');
addpath('../MDSadata');
% adtest(); 
addpath('../tests');
test01();
%Grat_sigma_temporal(); 
 
end 
function verybasictest()
% load('ptest');
if nargin < 1, 
    M = 3; 
    C = zeros(M,M);
    C(2,1) = 0.3;  
    C(1,2) = -C(2,1);       
    C(3,1) = .9;
    
    C = C - eye(M)/4;
    dT = 0.2; 
    
%     D = 
    T = 240;
    v = ones(1,T);     
    
    p.C = C;
    p.Pw.D = eye(M)/10^2;
    p.Pw.mu = zeros(M,1); % fix at zero.
    
    p.Pe.D = eye(M)/20^2;
    p.Pe.mu = zeros(M,1); % fix at zero for now. 
         
    p.dT = dT;
    p.s(1).v = v; 
    
    Phi_mu = spm_hrf(p.dT);
    
    pp = 3;  
    p.s(1).Phi = repmat(Phi_mu',[pp,1]);
%     p.s(1).Phi = p.s(1).Phi*0;
%     KK = 4;
%     p.s(1).Phi(:,1:KK) = 1/KK; 
    p.s(1).B = eye(M,pp);    
    
    p.s(1).u = zeros(M,T);
    p.s(1).U = zeros(M);
    
end 
% end
[p,y] = MDS_adata(p);

p.s(1).dudt = [p.s(1).u(:,2:end)-p.s(1).u(:,1:end-1),zeros(M,1)];

C0 = p.C;
TT = 100;
CC = zeros(M*M,TT);
for t=1:TT,
    p = MDS_sample_C(p);
    p = MDS_sample_s(p);
    %%
    subplot(4,1,4);
    %imagesc( [p.C,C0]); colorbar; 
    %pause(0.01);
    CC(:,t) = p.C(:);
    hold off;
    plot(CC(:,1:t)'); hold all;
%     tt = 1:t;
    plot(t, (mean(CC(:,1:t),2) ), 'o' );
    plot(t, C0(:),'.');
    hold off;
    II = [min(C0(:)),max(C0(:))];
    
    plot(II,II,'k-');  
    hold on;
    plot( C0(:), (mean(CC(:,1:t),2) ), 'o' ); axis equal;
    
    
    
    
    
    drawnow; 
end


 end
