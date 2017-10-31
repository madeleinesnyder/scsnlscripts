function p = MDS_fward(p,sub)
if nargin < 1, test(); return; end
if nargin < 2, 
    S = 1:length(p.s); 
else 
    S = sub;
end
for k=S,
    p = mds_sub_forward(p,k); % consider returning only the sub.
end
end

function p = mds_sub_forward(p,sub)
p.s(sub).x = s2x(p, p.s(sub).s);

%%
p.s(sub).yhat = permute( sum(bsxfun(@times, p.s(sub).B', mmx('mult',p.s(sub).Phi,p.s(sub).x)),1),[2,3,1]  ); 
%%
p.s(sub).xf = s2x(p, p.s(sub).sf); 

p.s(sub).yhatf = permute( sum(bsxfun(@times, p.s(sub).B', mmx('mult',p.s(sub).Phi,p.s(sub).xf)),1),[2,3,1]  ); 
if isfield(p.s(sub),'y'),
    p.s(sub).e = p.s(sub).y-p.s(sub).yhat-p.s(sub).yhatf;
end
%%  
M = size(p.C,1);
T = size(p.s(sub).y,2); 

s_tm1 =  [zeros(M,1), p.s(sub).s(:,1:end-1)];
if p.U_meanstim,
    s_tm1 = s_tm1 -  sum(bsxfun(@times, p.s(sub).U,  p.s(sub).u),3);
end

m = mmx('mult', p.dT * p.C, s_tm1 );
m = permute(m,[1,3,2]); 
v = permute(p.s(sub).v,[1,3,2])*1;
shat2 = permute(mmx('mult',m,v),[1,3,2]);
shat2(:,2:end) = shat2(:,2:end) + p.s(sub).s(:,1:end-1);
if ~p.U_meanstim,
    shat2 = shat2 +  sum(bsxfun(@times, p.s(sub).U, p.s(sub).u),3);
end
dw = p.s(sub).s-shat2;
p.s(sub).w = dw;  
%%  
if p.debug, 
    [S_inv_chol,q] = MDS_mk_S_matrix(p,sub);
    S_inv = S_inv_chol' * S_inv_chol;    
    lp1 = limvnpdf(vec(p.s(sub).s')',S_inv\q,S_inv);%+...
    if p.Pw_sigma_temporal,
        lp2 = 0; 
        for m=1:M,        
            lp2 = lp2 + limvnpdf(dw(m,:),zeros(T,1), p.Pw(m).S_inv );        
        end
    else
        lp2 = limvnpdf(dw',zeros(M,1), p.Pw.S_inv );    
    end
    
    if ~eqc(lp1,lp2)
       [lp1,lp2]
       assert(false);
    end
    %%
    assert(eqc(lp1,lp2) );
end
%%

if false, % who knows what this is. 
    %%
    dw = p.s(sub).s;
    dw(:,2:end) = diff(p.s(sub).s,1,2);
    for j=1:size(p.C,3),
        dw(:,2:end) = dw(:,2:end) - p.dT*p.C(:,:,j) * bsxfun(@times,p.s(sub).s(:,1:end-1),p.s(sub).v(1:end-1));
    end
    assert(eqc(p.s(sub).w,dw));    
end
end