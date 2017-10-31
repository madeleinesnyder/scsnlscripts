function x = s2x(p,s)
%%
L = size(p.s(1).Phi,2);
[K,T] = size(s);
x = zeros(L,K,T);

for t=1:T,
    dd = min(L-1,t-1);
    x(1:dd+1,:,t) = s(:,t:-1:t-dd)';
end

end