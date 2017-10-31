function sample_C_testpad()
clc;
M = 3; T = 10;
A = rand(M,T);
B = rand(M,T);
C = rand(M);
Se_inv = rand(M); Se_inv = Se_inv * Se_inv';

s0 = -1/2 * trace( (A-C*B)' * Se_inv * (A-C*B) )
s1 = -1/2 * trace(A' * Se_inv * A) + trace(A' * Se_inv * C * B)  ...
    -1/2 * trace(B' * C' * Se_inv * C * B)

Cv = C';
Cv = Cv(:);
SS_inv = kron(Se_inv, B * B');

s1 = -1/2 * trace(A' * Se_inv * A) + trace(B' * C' * Se_inv * A)  ...
    -1/2 * Cv' * SS_inv * Cv; %trace(B' * C' * S_inv * C * B)

q1 = kron(Se_inv * A,  eye(M)) * B(:);
m = B * A' * Se_inv;
q = m(:); % same_same;

s2 = -1/2 * trace(A' * Se_inv * A) + m' * Cv    ...
    -1/2 * Cv' * SS_inv * Cv; %trace(B' * C' * S_inv * C * B)


[s0,s1,s2] - s0
%%
end
