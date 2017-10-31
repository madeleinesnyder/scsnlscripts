function p = sample_s(p)
if nargin < 1, MDS(); return; end
% p = mds_fward(p);

if p.Pw.type == 2 && ~isfield(p.Pw, 'c'),
    T = size(p.s(1).w,2);
    for t=1:T
        [~, p.Pw.c(t).SW_inv, p.Pw.c(t).muProj] = GPmarg(p.Pw, p.s(1).w,t:min(T,t+1));
        p.Pw.c(t).SW_inv = inv(p.Pw.c(t).SW_inv);
        p.Pw.c(t).not_t1 = true(1,T);
        p.Pw.c(t).not_t1(t:min(T,t+1)) = false; 
    end
end


for k=1:length(p.s),
    if p.Pw.type == 1, 
        p.s(k) = dsample_s(p,k);
    else
        p.s(k) = dsample_s_GP(p,k);        
    end    
end 
if p.debug, MDS_selfcheck(p); end 
%% 
% p2 = mds_fward(p)
% p2.s.yhat - p.s.yhat
% p2.s.e - p.s.e

end

function ps = dsample_s_GP(p,k)      
C = p.C;
ps = p.s(k);
DLM = ps.B * ps.Phi; 
[M,T] = size(ps.y);
GPM = p.Pw.type == 2;
if p.debug,
    MDS_selfcheck(p);
end

GPM = p.Pw.type == 2; 
p.NWe.S_inv = p.NWe.S_inv;
% [~,SW,muProj] = GPmarg(p.Pw,ps.w,1:2);

    
%%
for t=1:T,
    if GPM,
        muProj = p.Pw.c(t).muProj;
        SW_inv =  p.Pw.c(t).SW_inv;
    end
    
    if p.debug,
        po = p;
        po.s(k) = ps;
        MDS_selfcheck(po);        
    end
    Ct = sum(bsxfun(@times, C, permute(ps.v(:,t), [2, 3, 1])),3);
    if t == 1, 
        st0 = ps.w(:,t)*0; 
    else
        st0 = Ct * ps.s(:,t-1);
    end    
    st0 = st0 + p.U * ps.v(:,t);    
    if t < T,
        Ct_1 = sum(bsxfun(@times, C, permute(ps.v(:,t+1), [2, 3, 1])),3);
        st_11 = ps.s(:,t+1) - Ct_1 * st0  - p.U * ps.v(:,t+1); 
        r = 1; 
    else
        st_11 = st0*0;
        Ct_1 = Ct*0; 
        r = 0; 
    end
    if GPM,
        MUW = ps.w(:,p.Pw.c(t).not_t1) * muProj; 
        if p.debug,
            [MUW2,SW2,~] = GPmarg(p.Pw,ps.w,t:t+r);            
            d_ = mean(mean(abs(MUW2 - MUW)));
            assert(d_ <10^-6);
            assert(eqc(inv(SW2),SW_inv));
            %%
        end
         
        mu_w = MUW;
%         SW_inv = inv(SW); %eye(size(MUW,1))*1/SW;
    else
        mu_w = p.Pw.mu;
        Sw_inv = p.Pw.S_inv;
    end
    psL = size(ps.Phi,2); 
    if ~GPM,
        st_11 = st_11 - r * mu_w(:,end);          
    end
    L = psL + min(0,T+1-t-psL);
    ll = 1:L;   
    DLMll = DLM(:,ll);
    Y = ps.y(:,t+ll-1) - ps.yhat(:,t+ll-1) + bsxfun(@times, DLMll, ps.s(:,t)-st0);
    Y = bsxfun(@minus, Y, p.NWe.mu);
    %%
%     clc
    Sc_inv = zeros(M); 
    K = 0;
    b = zeros(M,1);
    s0 = 0; 
    s1 = 0; 
    for i=1:M, % GMM
        a = st_11(i);
        A = zeros(2,1); 
        A(2) = a;
        A = A - MUW(i,:)';
        
        B = zeros(2,M);
        B(1,i) = -1; 
        B(2,:) = Ct_1(i,:);
        K = K + A' * SW_inv * A;
        Sc_inv = Sc_inv + B' * SW_inv * B;
        b = b + B' * SW_inv * A;
        
        wit = ps.w(i,t:t+r)'; 
        wt = ps.w(:,t);
        s0 = s0 + -1/2 * (wit - MUW(i,:)')' * SW_inv * (wit - MUW(i,:)');
        
    end
%     dde = bsxfun(@minus,ps.e(:,t+ll-1), p.NWe.mu);  
    
    %
%     s0 - (-1/2 * wit' * SW_inv * wit + B' * wt - 1/2 * K) %+ 1/2 * (A - B*wt)' * SW_inv * (A-B * wt)
%     s0 - (-1/2 * wt' * Sc_inv * wt + b' * wt - 1/2 * K)
    
    
    %%
%     s0 = s0 - 1/2 * trace(dde' *  p.NWe.S_inv * dde);
    %
    m = sum(DLMll .* (p.NWe.S_inv * Y),2);
    if ~GPM,
        Sc_inv = Sw_inv + r * Ct_1' * Sw_inv * Ct_1 + (DLMll * DLMll') .* p.NWe.S_inv; % NW
        b = (r*Ct_1'  * Sw_inv * st_11+ m) + Sw_inv * mu_w ;   % NW
    else        
        Sc_inv = Sc_inv + (DLMll * DLMll') .* p.NWe.S_inv; % GMM
        b = b + m; % + SW_inv * mu_w ;    % GMM
        K = K + trace(Y' * p.NWe.S_inv * Y);
    end
    
    
%     s0 - (-1/2 * wt' * Sc_inv * wt + b' * wt - 1/2 * K)
    %% independent test of the other factorization. 
    
    
    
    %%
    if ps.debug && false,        
        %% do debug code here.
        if GPM,
            ddw = bsxfun(@minus, ps.w, p.Pw.mu');
        else
            ddw = bsxfun(@minus, ps.w, p.Pw.mu);
        end
        dde = bsxfun(@minus, ps.e, p.NWe.mu);
        
        s0(t) = -1/2 *trace(ddw' * Sw_inv * ddw ) - 1/2 * trace(dde' *  p.NWe.S_inv * dde);    
        de = bsxfun(@minus,ps.e(:,setdiff(1:T, t+ll-1)), p.NWe.mu);   
        if GPM,
            dw =  bsxfun(@minus, ps.w(:,setdiff(1:T, t:t+1)), p.Pw.mu');                    
        else
            dw =  bsxfun(@minus, ps.w(:,setdiff(1:T, t:t+1)), p.Pw.mu);        
        end
        
        L0 =  -1/2 *trace(dw' * p.Pw.S_inv * dw) - 1/2 * trace(de' *  p.NWe.S_inv * de);
        L1 = L0 - 1/2 * r*st_11' * p.Pw.S_inv * st_11 - 1/2 * trace(Y' * p.NWe.S_inv * Y);        
        wt = ps.w(:,t);
        L1 = L1 - 1/2 * p.Pw.mu' * p.Pw.S_inv * p.Pw.mu;% - 1/2* 0 * p.NWe.mu' * p.NWe.S_inv * p.NWe.mu;
        
        s1(t) = L1 - 1/2 * (wt - Sc_inv\b)' * Sc_inv * (wt - Sc_inv\b) + 1/2 * b' * (Sc_inv\b);            
        if abs(s0(t) -s1(t)) / abs(s0(t)+s1(t)) > 10^-10,    
            assert(false);
        end  
    end 
    w_new = imvnrnd(Sc_inv\b, Sc_inv)';
    %% update cache. 
    ps.yhat(:,t+ll-1) = ps.yhat(:,t+ll-1) +  bsxfun(@times, DLMll,w_new-ps.w(:,t));
    
    ps.s(:,t) = ps.s(:,t) - ps.w(:,t) + w_new;
    ps.w(:,t) = w_new; 
    if t < T,
        ps.w(:,t+1) = ps.s(:,t+1)-Ct_1 * ps.s(:,t) - p.U * ps.v(:,t+1);
    end
    ps.e = ps.y - ps.yhat;         
    for i=ll,
        ps.x(i, :, t+i-1) = ps.s(:,t);
    end
    if p.debug,
        x1 = ps.x(:, :, t+ll-1);
        x2 = mtimesx(ps.s2x(:,:,t+ll-1),ps.s,'t');        
        assert(sum(abs(x1(:)-x2(:))) < 10^-10);
    end     
    if p.debug, 
        p2 = p; 
        p2.s(k) = ps; 
        MDS_selfcheck(p2);
        %%
        MDS_dlp_test(p2,po,w_new, po.s(k).w(:,t),Sc_inv\b,Sc_inv);    
    end
end 
end

function ps = dsample_s(p,k)      
C = p.C;
ps = p.s(k);
DLM = ps.B * ps.Phi; 
[~,T] = size(ps.y);
GPM = p.Pw.type == 2;
if p.debug,
MDS_selfcheck(p);
end
%%
for t=1:T,
    if p.debug,
        po = p;
        po.s(k) = ps;
        MDS_selfcheck(po);        
    end
    Ct = sum(bsxfun(@times, C, permute(ps.v(:,t), [2, 3, 1])),3);
    if t == 1, 
        st0 = ps.w(:,t)*0; 
    else
        st0 = Ct * ps.s(:,t-1);
    end    
    st0 = st0 + p.U * ps.v(:,t);    
    if t < T,
        Ct_1 = sum(bsxfun(@times, C, permute(ps.v(:,t+1), [2, 3, 1])),3);
        st_11 = ps.s(:,t+1) - Ct_1 * st0  - p.U * ps.v(:,t+1); 
        r = 1; 
    else
        st_11 = st0*0;
        Ct_1 = Ct*0; 
        r = 0; 
    end
    if GPM,
        [MUW,SW] = GPmarg(p.Pw,ps.w,t);
        mu_w = MUW;
        Sw_inv = eye(size(MUW,1))*1/SW;
    else
        mu_w = p.Pw.mu;
        Sw_inv = p.Pw.S_inv;
    end
    psL = size(ps.Phi,2); 
    st_11 = st_11 - r * mu_w;     
     
    L = psL + min(0,T+1-t-psL); 
    ll = 1:L;   
    DLMll = DLM(:,ll);
    Y = ps.y(:,t+ll-1) - ps.yhat(:,t+ll-1) + bsxfun(@times, DLMll, ps.s(:,t)-st0);
    Y = bsxfun(@minus, Y, p.NWe.mu);
    Sc_inv = Sw_inv + r * Ct_1' * Sw_inv * Ct_1 + (DLMll * DLMll') .* p.NWe.S_inv; %+ diag(sum(DLMll.* (p.NWe.S_inv * DLMll),2));
    m = sum(DLMll .* (p.NWe.S_inv * Y),2);
    b = (r*Ct_1'  * Sw_inv * st_11+ m) + Sw_inv * mu_w ;    
    %%
    if ps.debug && false,        
        %% do debug code here.
        if GPM,
            ddw = bsxfun(@minus, ps.w, p.Pw.mu');
        else
            ddw = bsxfun(@minus, ps.w, p.Pw.mu);
        end
        dde = bsxfun(@minus, ps.e, p.NWe.mu);
        
        s0(t) = -1/2 *trace(ddw' * Sw_inv * ddw ) - 1/2 * trace(dde' *  p.NWe.S_inv * dde);    
        de = bsxfun(@minus,ps.e(:,setdiff(1:T, t+ll-1)), p.NWe.mu);   
        if GPM,
            dw =  bsxfun(@minus, ps.w(:,setdiff(1:T, t:t+1)), p.Pw.mu');                    
        else
            dw =  bsxfun(@minus, ps.w(:,setdiff(1:T, t:t+1)), p.Pw.mu);        
        end
        
        L0 =  -1/2 *trace(dw' * p.Pw.S_inv * dw) - 1/2 * trace(de' *  p.NWe.S_inv * de);
        L1 = L0 - 1/2 * r*st_11' * p.Pw.S_inv * st_11 - 1/2 * trace(Y' * p.NWe.S_inv * Y);        
        wt = ps.w(:,t);
        L1 = L1 - 1/2 * p.Pw.mu' * p.Pw.S_inv * p.Pw.mu;% - 1/2* 0 * p.NWe.mu' * p.NWe.S_inv * p.NWe.mu;
        
        s1(t) = L1 - 1/2 * (wt - Sc_inv\b)' * Sc_inv * (wt - Sc_inv\b) + 1/2 * b' * (Sc_inv\b);            
        if abs(s0(t) -s1(t)) / abs(s0(t)+s1(t)) > 10^-10,    
            assert(false);
        end  
    end 
    w_new = imvnrnd(Sc_inv\b, Sc_inv)';
%     w_new = Sc_inv\b;%, Sc_inv)';
    
    
    %% update cache. 
    ps.yhat(:,t+ll-1) = ps.yhat(:,t+ll-1) +  bsxfun(@times, DLMll,w_new-ps.w(:,t));
    
    ps.s(:,t) = ps.s(:,t) - ps.w(:,t) + w_new;
    ps.w(:,t) = w_new; 
    if t < T,
        ps.w(:,t+1) = ps.s(:,t+1)-Ct_1 * ps.s(:,t) - p.U * ps.v(:,t+1);
    end
    ps.e = ps.y - ps.yhat;         
    for i=ll,
        ps.x(i, :, t+i-1) = ps.s(:,t);
    end
    if p.debug,
        x1 = ps.x(:, :, t+ll-1);
        x2 = mmx('mult',ps.s2x(:,:,t+ll-1),ps.s,'nt');        
        assert(sum(abs(x1(:)-x2(:))) < 10^-10);
    end     
    if p.debug,
        p2 = p; 
        p2.s(k) = ps; 
        MDS_selfcheck(p2);
        MDS_dlp_test(p2,po,w_new, po.s(k).w(:,t),Sc_inv\b,Sc_inv);    
    end
end 
end

% work out alternative prior structure. 
