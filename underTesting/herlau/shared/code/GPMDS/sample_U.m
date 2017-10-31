function p = sample_U(p,opts)
if nargin < 1, test(); end
% sample the COLUMNS of U.
J = size(p.C,3);
for j=1:J,    
    p = dsample_U(p,j);
end

end
function p = dsample_U(p,i)
4;
%%  
S_inv = p.NWU.S_inv;
q = p.NWU.S_inv * p.NWU.mu; 
Ui = p.U(:,i);
if p.debug,     
    dU = Ui - p.NWU.mu;
    s0 = -1/2 * dU' * p.NWU.S_inv * dU; 
    K = p.NWU.mu' * p.NWU.S_inv * p.NWU.mu;     
    po = p;   
end
%
for k=1:length(p.s),
    sh = bsxfun(@minus, p.s(k).w, p.NWw.mu) + p.U(:,i) * p.s(k).v(i,:);       
    S_inv = S_inv + p.NWw.S_inv * sum(p.s(k).v(i,:).^2);
    q = q + p.NWw.S_inv * (sh * p.s(k).v(i,:)');

    if p.debug,
        dw = bsxfun(@minus, p.s(k).w,p.NWw.mu);
        s0 = s0 -1/2 * trace(dw' * p.NWw.S_inv * dw);
        K = K  +  trace(sh' * p.NWw.S_inv * sh);
    end
end 
U_new = imvnrnd(S_inv\q, S_inv)';
if p.debug,
    s1 = -1/2 * Ui' * S_inv * Ui + q' * Ui - 1/2 * K; % 1/2 * q' * (S_inv\q) - 1/2*K
    if abs(s1-s0) / abs(s1+s0) > 10^-10,
        dw1 = bsxfun(@minus, p.s(1).w, p.NWw.mu);
        dw2 = bsxfun(@minus, sh, Ui * p.s(k).v(i,:));
        disp(dw1-dw2)
        assert(false);
    end
end
p.U(:,i) = U_new;
% reset other constants.  
for k=1:length(p.s)
    p.s(k).w = p.s(k).w + (-U_new+Ui) * p.s(k).v(i,:);
%     p.s(k).w = bsxfun(@plus, p.s(k).w, p.NWw.mu);
end
% p2 = MDS_fward(p);
% p.s(k).w-p2.s(k).w
%% 
if p.debug,
    MDS_selfcheck(p);
    MDS_dlp_test(p,po, U_new, po.U(:,i), S_inv\q, S_inv);    
end

%%
end