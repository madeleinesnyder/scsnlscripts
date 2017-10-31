function s2x = mk_s2x(T,L)
% s2x = NaN; return;
s2x = [-1, 1];
return;
s2x = zeros(L,T,T);
for t=1:T,
    m0 = eye(L); m0 = m0(:,end:-1:1);
    m0 = [m0 ; zeros(max(0,T-L),L)]';
    m0 = m0 * lagm(T,L-t);
    s2x(:,:,t) = m0; 
end
    

end