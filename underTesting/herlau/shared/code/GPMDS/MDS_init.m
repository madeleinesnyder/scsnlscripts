function p = MDS_init(p,opts)
if nargin < 1, test(); return; end
cdir = fileparts(mfilename('fullpath'));
addpath(fullfile(cdir, './dists'));

cdir = fileparts(mfilename('fullpath')); 
addpath(fullfile(cdir, 'mmx'))  
mmx(1);
if nargin < 1, test(); return; end
pdef = struct(); 
pdef.debug = false; 

odef.T = 200; odef.S = 1; 

if nargin < 2, opts = struct(); end
if ~isfield(p,'dT'), p.dT = 0.75; end

%odef.GP_hrf = true;
% odef.Pe.kappa = 1;
odef.PPhi.Tscale = 10 / p.dT; % in seconds.   
odef.PPhi.sigf = 0.05 * p.dT; 
odef.PPhi.sign = 0.001 * p.dT;
odef.PPhi.type = 2; % GP. 

%%
odef.PU.type = 7; 
odef.PU.alpha = 1; 
odef.PU.beta = 1; 
odef.PU.opts = struct();
odef.PU.U_meanstim = false; 
odef.use_du = false;

%% P(e)
odef.Pe.type = 2;
odef.Pe.gptype = 2; 
odef.Pe.sigf = 1/3;
odef.Pe.sign = 1/20;
odef.Pe.l = 2; 
odef.Pe.opts = struct();

%% P(w) this does not seem super hot. 
odef.Pw.type = 5;
odef.Pw.sigma = 1/3;
odef.Pw.nu = 5;
odef.Pw.A = 1; 
odef.Pw.alpha = 1; 
odef.Pw.beta = 5; 
odef.Pw.opts = struct(); 
odef.Pw.opts.Tsigma = 4; 
odef.Pw.sigma_temporal = false; 
 
%% P(C)
odef.PC.type = 7; 
odef.PC.opts = struct();
odef.PC.alpha = 1; 
odef.PC.beta = 10; 
%% P(s)
odef.Psf.opts.Psign.alpha = 1; 
odef.PSf.opts.Psign.beta = 2;
odef.Psf.opts.Tl = 0; 
odef.Psf.opts.Tsigf = 0; 
odef.Psf.l = 1; 
odef.Psf.sign = 1; 
odef.Psf.sigf = 0; 
odef.Psf.gptype = 1; 
odef.Psf.USE_Psf = false;
%%
if isfield(p,'s'),
    odef.S = length(p.s);
    if isfield(p.s(1),'v'),odef.T = size(p.s(1).v,2); end
    if isfield(p.s(1), 'y'),
        odef.T = size(p.s(1).y,2);
    end
end 
%%
opts = ssfr(odef,opts);
GP_hrf = opts.PPhi.type == 2; 

[M,~,J] = size(p.C);
if ~isfield(opts,'p'), opts.p = M; end

% close all;
if opts.PPhi.type == 2, % GP. 
    Phi_mu = spm_hrf(p.dT); 
    L = length(Phi_mu);
end
%%
pdef.C = repmat(eye(M)/2, [1,1,J]) + rand(M,M,J)/100;
for k=1:length(p.s),
    T = size(p.s(k).v,2); %
    pdef.s(k).u = zeros(M,T);
    H = 1; 
    if isfield(p.s(k), 'u'), H = size(p.s(k).u,3); end
    pdef.s(k).U = zeros(M,1,H);
    if ~GP_hrf, 
        [Phi] = get_model_hrf_3basis(p.dT);
        pdef.s(k).Phi = Phi;
    end
end
p = ssfr(pdef,p);     
pss = p.s; 
p = rmfield(p,'s'); 
pdef = rmfield(pdef,'s');
if isfield(pss(1),'Phi') && ~GP_hrf, [pp,L] = size(pss(1).Phi); end

pdef.PB = mk_HT(1,true);
for h=1:size(pss(1).u,3),
    pdef.PU(h) = th_NGAM(zeros(M,1), opts.PU.alpha,opts.PU.beta,opts.PU.opts);
end

if opts.Pw.sigma_temporal,
    for m=1:M,
    pdef.Pw(m) = th_NGAM(zeros(T,1), opts.Pw.alpha, opts.Pw.beta, opts.Pw.opts); % you want a lot of obs. 
    end
else
    pdef.Pw = th_NGAM(zeros(M,1), opts.Pw.alpha, opts.Pw.beta, opts.Pw.opts);    
end
pdef.Pw_sigma_temporal = opts.Pw.sigma_temporal;
pdef.U_meanstim = opts.PU.U_meanstim; 

if opts.Pe.type == 5,
    pdef.Pe = mk_HT(M,false,opts.Pe.sigma);
else
    for k=1:M,
        %% P(e)
        pdef.Pe(k) = th_GP(opts.Pe.gptype,zeros(T,1),opts.Pe.l, opts.Pe.sigf, opts.Pe.sign, opts.Pe.opts);
        pdef.Pe(k).mu(:) = 0;        
        %% P(sf)
        pdef.Psf(k) = th_GP(opts.Psf.gptype,zeros(T,1),opts.Psf.l, opts.Psf.sigf, opts.Psf.sign, opts.Psf.opts);
    end
    p.USE_Psf = opts.Psf.USE_Psf;    
end

if opts.PC(1).type == 1, 
    pdef.PC = mk_HT(M*M,true);
end
 for j=1:J,
    if opts.PC(1).type == 6,
        pdef.PC = th_DL(zeros(M*M*J,1), opts.PC.a);
        pdef.PC = pdef.PC.MCMC(vec(p.C)' );
    elseif opts.PC(1).type == 7,   
        %% 
        jj = 1; 
        if length(opts.PC) > 1, jj = j; end        
        pdef.PC(j) = th_NGAM(zeros(M*M,1), opts.PC(jj).alpha, opts.PC(jj).beta,opts.PC(jj).opts );
        pdef.PC(j) = pdef.PC(j).MCMC(vec(p.C(:,:,j))',1);
     
    end
end  

if false,
    if opts.PC.type == 1,
        pdef.PC = mk_NW(M); 
    elseif opts.PC.type == 3,
         pdef.PC = mk_GM([M,M,J]);
    else
        assert(false);
    end
end

p = ssfr(pdef,p);
%%

if GP_hrf,
    p.PPhi = mk_GP(opts.p, L, opts.PPhi.Tscale, opts.PPhi.sign, opts.PPhi.sigf, true,Phi_mu);
    %% 
    for k=1:length(pss),
        pss(k).Phi = GP_rnd(p.PPhi)';
        pss(k).Phi = repmat(Phi_mu',[size(pss(k).Phi,1), 1]);
        pss(k).B = eye(M,M);        
    end
else
    for k=1:length(pss),
        pss(k).B = ones(M,pp)/10;
        pss(k).B(:,1) = 1;
    end
end
for sub=1:length(pss) 
    sd = pss(sub); 
    
    if opts.use_du && ~p.U_meanstim,
        sd.u = [zeros(size(sd.u,1),1,size(sd.u,3)), sd.u(:,2:end,:)-sd.u(:,1:end-1,:)];
    end    
    sd.s2x = mk_s2x(T,L);
    
    if ~isfield(sd,'y'),        
    end
    if isfield(sd,'y'),
        T = size(sd.y,2);
    end
    p2 = p;
    if isfield(p2, 's'), p2 = rmfield(p2,'s'); end
    p2.s(1) = sd;
    [S_inv_chol,q] = MDS_mk_S_matrix(p2,1);
    S_inv = S_inv_chol' * S_inv_chol;
    mu = S_inv\q; 
         
    if ~isfield(sd,'s'),
        s = imvnrnd(mu,S_inv);  
        s = reshape(s,[T,M])';
        sd.s = s; 
    end
    sd.sf = 0*rand(M,T)/100+zeros(M,T);
     
    sd.shat = reshape(mu,[T,M])';     

    if ~isfield(sd,'y'),
        x = mmx('mult',sd.s2x, sd.s, 'nt');
        sd.yhat = permute( sum(bsxfun(@times, sd.B', mmx('mult',sd.Phi,x)),1),[2,3,1]  );          
        for k=1:length(p.Pe),
            sd.y(k,:) = sd.yhat(k,:) + p.Pe(k).sample; %mvnrnd(p.Pe(m).mu,p.Pe(m).S,T)';        
        end        
    end
    %%
    p.s(sub) = ssfr(sd,pss(sub));            
end
end