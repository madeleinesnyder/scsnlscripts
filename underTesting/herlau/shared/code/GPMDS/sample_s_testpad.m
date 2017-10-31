function sample_s_testpad(p),
% derive the full marginal. 
rng(1);
clc;
T = 2; 
M = 3; 
J = 2;   
% mangler flere ting: Normalisering, HRF, skift i Ct. 
dT = 0.05; 
dT = 0.02; 
dT = 0.25; 

T = round(T / dT); 
disp(T);
D = eye(M) * 0.01; 
D = rand(M); D = D*D';
C = rand(M,M,J)*.5;
C = bsxfun(@minus, C, eye(M));

v = rand(J,T);
% v = v * 0 + 1; 


if J == 1, v = v*0+1; end
if false,
    for j=1:J,
        C = C * 0;
        dd = .9; 
        for k=1:M,
            C(k,k,j) = -0.2; 
            C(mod(k,M)+1,k,j) = dd; 
        end
        C(1,end,j) = -C(1,end);     
    end    
end 
% C = C + rand(M,M,J)/100;
Ct = bsxfun(@plus, C * dT, eye(M));
Ct2 = Ct;
Ct = zeros(M,M,T-1);
for t=1:T-1,
    for j=1:J,
        Ct(:,:,t) = Ct(:,:,t) + Ct2(:,:,j)*v(j,t);
    end
end
Ct = cat(3,zeros(M),Ct);


%SD_inv = rand(M,M); SD_inv = SD_inv * SD_inv'  + eye(M) /10;
SD_inv = inv(D * dT);
s = rand(M,T);
GS_inv = rand(T); GS_inv = GS_inv  * GS_inv' + eye(T)/10; 

% dT = 1;
Tscale = pi; % scale over which this one is smooth. 
sign = 1;  % turn off GP.
sigf = 1;   
GP = mk_GP(M,T,dT,Tscale, sign, sigf);
GS_inv = GP.S_inv; %inv( K(GP,(1:T)') );
s0 = zeros(M,1);
y = zeros(M,T); 
% d = 0.8;
% Ct(1,2) = d;
% Ct(2,3) = d;
% Ct(3,1) = d;
% lambda = 1;
% Ct(2,1) = -0.6;
% for t=1:T,
%     
% end 

sz0 = -1/2 * trace(s * GS_inv * s'); 
s0 = sz0;

dd0 = [];
for t=1:T,
    if t == 1, dd0(:,t) = s(:,t);
    else dd0(:,t) = s(:,t) - Ct(:,:,t) * s(:,t-1); end
    s0 = s0 -1/2 * trace( dd0(:,t)' * SD_inv * dd0(:,t));
end

Rf = zeros(T-1,T-1);
Rb = zeros(T,T-1);
for t=1:T-1,
    for j=t:T,
        Rf(t+1,t) = 1; 
        Rb(t,t+1) = 1; 
    end
end 
Rf = eye(T);

if J == 1 && false, 
    dd = s * Rf - Ct(:,:,1) * s * Rb;
    s1 = -1/2 * trace( dd' * SD_inv * dd);

    s2 = -1/2 * trace( (Rf' * s' -Rb' * s' * Ct') * SD_inv * (s * Rf - Ct * s * Rb));

    s3 = -1/2 * trace( (Rf' * s' ) * SD_inv * (s * Rf ));
    s3 = s3 + trace( Rf' * s' * SD_inv *  Ct * s * Rb);
    s3b = s3 + sz0; 

    s3 = s3 - 1/2 * trace( ( Rb' * s' * Ct') * SD_inv * ( Ct * s * Rb));
    H =  Ct' * SD_inv *  Ct;
    X = s';
    vX = X(:);

    s4 = s3b - 1/2 * trace(  Rb' * X * H * X' * Rb );


    % [s0, s1, s2,s3,s4] - s0


    s6 = sz0 -1/2 * vX' * kron( Ct' * SD_inv *  Ct, Rb * Rb') * vX ...
        - 1/2 * vX' * kron( SD_inv, Rf * Rf') * vX ...
        + vX' * kron(SD_inv *  Ct, Rf * Rb'  ) * vX;
    % [s0, s1, s2,s3,s4,s6] - s0
    s7 = -1/2 * vX' * kron(eye(M), GS_inv) * vX -1/2 * vX' * kron( Ct' * SD_inv *  Ct, Rb * Rb') * vX ...
        - 1/2 * vX' * kron( SD_inv, Rf * Rf') * vX ...
        + vX' * kron(SD_inv *  Ct, Rf * Rb'  ) * vX;

    SS_inv = kron(eye(M), GS_inv) + kron( Ct' * SD_inv *  Ct, Rb * Rb') ...
     + kron( SD_inv, Rf * Rf') - kron(SD_inv *  Ct, Rf * Rb'  ) - kron(Ct' * SD_inv,Rb * Rf' );

end

m0 = zeros(T);
if J > 1 || true, 
    SS_inv4 = kron(eye(M), GS_inv) ;
    SS_GP = kron(eye(M), GS_inv) ;
    SS_inv5 = SS_inv4; 
    SS_inv6 = SS_inv4; 
    SS_inv7 = SS_inv4;
    SS_inv8 = zeros(size(SS_inv4));
    SS_inv81 = zeros(size(SS_inv4));
    
    SS_inv9 =  zeros(size(SS_inv4));
    
    L = chol(SD_inv);
%     L' * L - SD_inv
%     SS_inv10 = SS_inv4; 
    for t=1:T, 
        et = ((1:T)' == t)*1;       
        
        p = M; q = M; 
        r = 1; ss = T;
        Spr = PS(p,r);
%         Spq = PS(p,q);
%         PS(4,4)
        Sqs = PS(q,ss);
        
        LT = L * Spr';
%         LT= LT + rand(M);
%         dd3 = kron(eye(T),LT) + kron(Rb, Ct(:,:,end));
        etp = (1:T)' == t-1;
        Ctt = Ct(:,:,end);
        
        dd8 = L *(kron(Ctt,etp' ) - kron(eye(M),et')) ;
        dd8b = LT * (kron(etp',Ct(:,:,end) ) - kron(et',eye(M)))* Sqs;  
        dd8c = LT * (kron(etp',Ct(:,:,end) ) - kron(et',eye(M)));
        
        dd8-dd8b;
        
        SS_inv81 = SS_inv81 + dd8b' * dd8b;
        SS_inv8 = SS_inv8 + dd8c' * dd8c;                
        
        m1 = dd8b' * dd8b;
        m2 = Sqs' * dd8c' * dd8c * Sqs;
        %%
        clc
        ll = kron(eye(T),LT) - kron( (etp * et')',LT * Ct(:,:,end));
        dd8c' * dd8c;
        ll' * ll;
        SS_inv9 = SS_inv9 + ll' * ll;
        %%        
        dd7 = kron(L, et')*(kron(Ctt,Rb') - kron(eye(M),eye(T)));
        dd7 = kron(L* Ctt, etp') - kron(L, et') * kron(eye(M),eye(T));
        dd7 = L * kron(Ctt, etp') - L * kron(eye(M),et');
        dd7 = L*(kron(Ctt, etp') -  kron(eye(M),et'));
        
        
        SS_inv7 = SS_inv7 + dd7' * dd7;
         
        dd = kron(Ct(:,:,end), (Rb )' ) -kron(eye(M), eye(T));
        SS_inv4 = SS_inv4 + dd' * kron(SD_inv, et * et') * dd;    
         
         
%         SS_inv5 = SS_inv5 + 
        
        m0 = m0 + Rb * et * (Rb * et)';
        SS_inv5 = SS_inv5 + kron( SD_inv, et * et') + kron( Ct(:,:,t)' * SD_inv * Ct(:,:,t), Rb * et * (Rb * et)')- 2*kron( SD_inv * Ct(:,:,t),  et * (Rb * et)' ); 
    end
    %% THIS IS SUFFICIENT!
    L11 = kron(eye(T),LT) - kron( Rb',LT * Ct(:,:,end));
    [I,~] = find(Sqs);
    L11 = L11(:,I);    
    SS_inv11 = SS_GP + L11' * L11;
        
    % REAL DEAL for J>1.
    L12 = kron(eye(T),LT);
    Ct12 = bsxfun(@plus, C * dT, eye(M));    
    for j=1:J,
        Rbv = Rb;
        Rbv(Rb ~= 0) = v(j,1:end-1);
        L12 = L12 - kron( Rbv',LT * Ct12(:,:,j));        
    end     
    [I,~] = find(Sqs);
    L12 = L12(:,I);    
    SS_inv12 = SS_GP + L12' * L12;        

    %%
    SS_inv8 = SS_GP +  Sqs' * SS_inv8 * Sqs;
    SS_inv81 = SS_GP +  SS_inv81;
    SS_inv9 = SS_GP + Sqs' * SS_inv9 * Sqs;
    
    dd = kron(Ct(:,:,end), Rb') - kron(eye(M),Rf');
    SS_inv6 = SS_inv6 + dd' *kron(SD_inv, eye(T)) * dd;

%     SS_inv8 =  kron(eye(M), GS_inv) +     Sqs' * dd3' * dd3 * Sqs;
        
    SS_inv7 - SS_inv6;
%     SS_inv8 - SS_inv7
    
    X = s';
    vX = X(:);
%     SS_inv8 
    
    s8b = -1/2 * vX' * SS_inv4 * vX;
    s8c = -1/2 * vX' * SS_inv5 * vX;
    s9 =  -1/2 * vX' * SS_inv6 * vX;
    s10 =  -1/2 * vX' * SS_inv7 * vX;
    s11 =  -1/2 * vX' * SS_inv8 * vX;
    s12 =  -1/2 * vX' * SS_inv81 * vX;
%     s13 =  -1/2 * vX' * SS_inv9 * vX;
    s14 =  -1/2 * vX' * SS_inv11 * vX;
    s15 =  -1/2 * vX' * SS_inv12 * vX;
    
    
%     [SS_inv4 - SS_inv6]
    [s0,s8b, s8c,s10,s11,s12,s14,s15]-s0    
    [s0,s15]-s0    
    
%     return;
end
clc; close all;
%% neeearly done.
tic();
IT = 1:T;
L = M+1;
H = rand(M,L);
Se_inv = diag(diag(rand(M)));
Se_inv = Se_inv * Se_inv';
y = rand(M,T);
muy = rand(M,1);
s0Y = 0;
for t=IT,
    sl = zeros(M,L);
    for l=1:min(t, L),
        sl(:,l) = s(:,t-l+1);
    end     
    dd = y(:,t) - sum(H .* sl,2) - muy;
    s0Y = s0Y - 1/2 * (dd' * Se_inv * dd);
end

y2 = bsxfun(@minus,y,muy);
S_invY = zeros(T*M,T*M);
q = zeros(T*M,1);
for k=1:M, 
    b = H(k,end:-1:1)';
    BB = (b * b') * Se_inv(k,k);
    for t=IT
        Lmin = min(t,L);
        I = (t-Lmin+1:t) + T*(k-1);
        S_invY( I,I) = S_invY(I,I) + BB(L-Lmin+ (1:Lmin),L-Lmin+(1:Lmin));
        q(I) = q(I) + Se_inv(k,k) * b(L-Lmin+ (1:Lmin)) * y2(k,t);
    end
end

X = s';
Xv = X(:);
s1 = -1/2*Xv' * S_invY * Xv - 1/2 * trace(y2' * Se_inv * y2) + q' * Xv;
[s0Y,s1]-s0Y
toc()
%% Joint testing:
ss0 = s0 + s0Y;
S_inv = S_invY + SS_inv12;
s1 = -1/2*Xv' * S_inv * Xv - 1/2 * trace(y2' * Se_inv * y2) + q' * Xv;

L = cholcov(S_inv);
dd = L * (Xv-S_inv\q);
dq = L'\q; 
s2 = -1/2*(Xv-S_inv\q)' * S_inv * (Xv-S_inv\q) - 1/2 * trace(y2' * Se_inv * y2) + 1/2 * q' * inv(S_inv) * q;
s3 = -1/2 * (dd' * dd)  - 1/2 * trace(y2' * Se_inv * y2) + 1/2 * (dq' * dq); 
[ss0,s1,s2,s3]-ss0
%%
clc
% do a little trace test.
A = rand(10,5); X = rand(10,6); B = rand(6,5);

trace(A' * X * B);
Av = A(:);
Xv = X(:);
XTv = X'; XTv = XTv(:);

trace(A' * X * B) - vec(B *A')' * vec(X')
trace(A' * X * B) - vec(A*B')' * vec(X)


end
function M = vec(M),
M = M(:);
end
function cholsmp(x,C,SG_inv,SD_inv)
Rf = zeros(T-1,T-1);
Rb = zeros(T,T-1);
for t=1:T-1,
    for j=t:T,
        Rf(t+1,t) = 1; 
        Rb(t,t) = 1; 
    end
end

end
function S = PS(m,n)
In = eye(n);

for i=1:m,
    ei = (1:m)' == i;
    dm =     kron(ei', kron(In,ei)); 
    if i == 1,
        S = dm;
    else
        S = S + dm;
    end

end
end

