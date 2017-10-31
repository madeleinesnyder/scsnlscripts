function lp = MDS_logp(p)
if nargin < 1, test(); return; end
lp = 0; 
%%
if ~p.Pw_sigma_temporal      
    lp = lp + p.Pw.logp([p.s.w]');    
end

S = length(p.s);
[M,T] = size(p.s(1).y);
for m=1:M,
    E = zeros(S,T); 
    W = zeros(S,T);
    SF = zeros(S,T);    
%     UU = zeros(S,T);
    for s=1:S,
        E(s,:) = p.s(s).e(m,:);
        W(s,:) = p.s(s).w(m,:);
        SF(s,:) = p.s(s).sf(m,:);               
    end  
    lp = lp + p.Pe(m).logp(E);
    if p.Pw_sigma_temporal,
        lp = lp + p.Pw(m).logp(W); 
    end    
    if p.USE_Psf,
        lp = lp + p.Psf(m).logp(SF);        
    end    
end 
%%
for k=1:length(p.s),     
    lp = lp + logData(p.PB,p.s(k).B(:) );    
    if isfield(p, 'PPhi'),
        lp = lp + logData(p.PPhi, p.s(k).Phi);        
    end  

end
for h=1:length(p.PU), % one for each relation. 
    U = zeros(S,M);
    for k=1:length(p.s), 
        U(k,:) = p.s(k).U(:,:,h);
    end
    lp = lp + p.PU(h).logp( U ); %logData(p.PU, p.s(k).U );
end 

[M,~,J] = size(p.C);
if isfield(p, 'PPhi'),    
    lp = lp + logPrior(p.PPhi);
end
lp = lp + logPrior(p.PB);
% lp = lp + logPrior(p.PU);
for j=1:J,
    lp = lp + p.PC(j).logp(vec(p.C(:,:,j))'); %1, 2, 3]))');
end

end