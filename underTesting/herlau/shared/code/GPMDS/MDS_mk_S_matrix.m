function [S_inv_chol,q,S_inv,mu,sUU] = MDS_mk_S_matrix(p,k,chol_decompose)
if nargin < 1, test(); return; end
if nargin < 3, chol_decompose = true; end
C = p.C;
[M,~,J] = size(C);
T = size(p.s(k).v,2);
U = p.s(k).U; 
v = p.s(k).v;
%% do you even need the shuffle matrices?
% if ~p.U_meanstim,
%     dudt = p.s(k).du;
% end
% Sw_inv = ; 
if ~isfield(p.s(1),'s'), p.debug = false; end

if p.debug,  % this seems kindof like a bad idea... shoudl just be p.s(k).w. 
    Ct = zeros(M,M,T-1);
    s = p.s(k).s;
    for t=1:T-1,
        Ct(:,:,t) = eye(M); 
        for ii=1:J,
            Ct(:,:,t) = Ct(:,:,t) + C(:,:,ii)*v(ii,t+1)*p.dT;
        end
    end
    Ct = cat(3,zeros(M),Ct);
    dd0 = []; 
    for t=1:T, 
        if ~p.U_meanstim, 
            UUSTIM = sum( bsxfun(@times,U, p.s(k).u(:,t,:) ), 3);
        else
            CCT = zeros(M);
            for j=1:J,
                CCT = CCT + p.C(:,:,j) * p.s(k).v(j,t);
            end
            UUSTIM = -p.dT * CCT * sum(bsxfun(@times, U, p.s(k).u(:,t,:) ),3);  
        end
        if t == 1, dd0(:,t) = s(:,t) - UUSTIM; 
        else dd0(:,t) = s(:,t) - Ct(:,:,t) * s(:,t-1) - UUSTIM; end        
    end 
    s0 = 0; 
    if p.Pw_sigma_temporal,
        for m=1:M,
            s0 = s0 + -1/2 * dd0(m,:) * p.Pw(m).S_inv * dd0(m,:)';    
        end
    else
        s0 = -1/2 * trace( dd0' * p.Pw.S_inv * dd0);        
    end 
    
end
Rb = lagm(T);

%% 
if p.U_meanstim,  
%     UU = zeros(M,T);
    J = size(p.C,3);        
    for j=1:J,
        V = diag(v(j,:));   
%         UU = UU -p.dT * p.C(:,:,j) * sum(bsxfun(@times, U, p.s(k).u),3) * V;            
    end
else
%     UU = sum(bsxfun(@times,U, p.s(k).u),3);
end
    

if p.Pw_sigma_temporal,
    A = cat(3, eye(M), -eye(M));
    B = cat(3, eye(T), Rb);
    UU = zeros(M,T);
    for j=1:J,            
        A(:,:,end+1) = -p.dT * p.C(:,:,j);    
        V = diag(v(j,:));        
        B(:,:,end+1) = Rb * V;     
        
        if p.U_meanstim, 
            UU = UU -p.dT * p.C(:,:,j) * sum(bsxfun(@times, U, p.s(k).u),3) * V;            
        end
    end    
    if ~p.U_meanstim,
        UU = sum(bsxfun(@times,U, p.s(k).u),3);
    end
    %%
    if p.debug,
        %%
        dww = sum( mmx('mult',mmx('mult',A, p.s(k).s),B),3) - UU;        
        assert(eqc( p.s(k).w , dww));         
    end
    
    %%
    S_inv = zeros(M*T);
    q = zeros(M*T,1);  
    sUU = 0; 
    for m=1:M, 
        L = cholcov(p.Pw(m).S_inv); 
        Wt = zeros(T,M*T);
        for j=1:size(A,3),                                      
            Wt = Wt + kron(A(m,:,j), L * B(:,:,j)' );                              
            q = q + vec( (UU(m,:) * p.Pw(m).S_inv * B(:,:,j)')' * A(m,:,j) );            
        end
        S_inv = S_inv + Wt' * Wt;       
        
        if nargout > 2 || p.debug
            sUU = sUU -1/2 * trace( UU(m,:) * p.Pw(m).S_inv *UU(m,:)' );
        end        
    end     
    %%    
    if chol_decompose || p.debug,
        S_inv_chol = cholcov(S_inv);    
    else
%        disp('Cholesky decomposition FAILED');
       S_inv_chol = 42+zeros(4,1);  
    end 
    if p.debug, 
        s1 = -1/2 * vec(s')' * S_inv_chol' * (S_inv_chol * vec(s') );            
        s1 = s1 + sUU; 
        s1 = s1 + q' * vec(s');
        %% recompute difference using A, B, UU.
        if ~eqc(s0,s1),
            ([s0,s1]-s0)/mean([s0,s1])
%             assert(false);
        end
    end
else  
    L = chol( p.Pw.S_inv );
    S_inv_chol = kron(L, eye(T))  - kron(L, Rb'); 
    UU = zeros(M,T);
    
    if p.U_meanstim,  
        J = size(p.C,3);        
        for j=1:J,
            V = diag(v(j,:));   
            UU = UU -p.dT * p.C(:,:,j) * sum(bsxfun(@times, U, p.s(k).u),3) * V;            
        end
    else
        UU = sum(bsxfun(@times,U, p.s(k).u),3);
    end
    
    q2 = UU' * p.Pw.S_inv  - Rb * UU' * p.Pw.S_inv ;        
    for j=1:J,
        V = diag(v(j,:));
        S_inv_chol = S_inv_chol - p.dT * kron(L * p.C(:,:,j), V' * Rb');
        q2 = q2 - p.dT * Rb * V * UU' * p.Pw.S_inv  * p.C(:,:,j);        
    end 
    q = vec(q2);  
    if (nargout > 2 || p.debug) 
        sUU = -1/2 * trace( (UU)' * p.Pw.S_inv  * (UU) );        
    end
    if any(~isfinite(S_inv_chol(:))),
       assert(false);
    end
end

%%l
if p.debug, 
    if ~p.Pw_sigma_temporal,
        vs = vec(s');
        S_inv_chol2 = kron(L, eye(T))  - kron(L, Rb'); 
%         UU = U * dudt;
        q2 = UU' * p.Pw.S_inv - Rb * UU' * p.Pw.S_inv;        
        for j=1:J,
            V = diag(v(j,:));
            S_inv_chol2 = S_inv_chol2 - p.dT * kron(L * p.C(:,:,j), V' * Rb');
            q2 = q2 - p.dT * Rb * V * UU' * p.Pw.S_inv * p.C(:,:,j);        
        end 
        q2 = vec(q2); 
        s3 = -1/2* vs' * (S_inv_chol2' * S_inv_chol2) * vs + q2' * vs - 1/2 * trace(UU' * p.Pw.S_inv * UU);

        s1q = q' *vs;
        SS_inv12 = S_inv_chol' * S_inv_chol;    
        s1 = -1/2 * vs' * SS_inv12 * vs + s1q + sUU;
    else
        s1 = -1/2 * (vec(s')' * S_inv_chol') * (S_inv_chol * vec(s') ) ;            
        s1 = s1 + sUU;
        s1 = s1 + q' * vec(s');
         
        s3 = s1;  
    end 
%     [s0, s1, s3] - s0 
    if ~eqc(s0,s1) || ~eqc(s0,s3),
         [s0, s1, s3] - s0
        assert(false);
%        34  
    end
end

if nargout > 2,
    if ~p.Pw_sigma_temporal,
        S_inv = S_inv_chol'*S_inv_chol;
    end    
end
if nargout > 3,
    mu = S_inv\q; 
end

end