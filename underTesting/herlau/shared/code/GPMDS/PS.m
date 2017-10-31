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