function R = lagm(n,l)
if nargin < 2, l = -1; end
R = ones(n);
R = tril(triu(R,-l),-l);
return;
ii = l+1:n;
jj = 1:n-l;
s = ii*0+1;
R2 = sparse(ii,jj,s,n,n);
% toc();
R = R2;
assert(all(R(:)==R2(:)));
end