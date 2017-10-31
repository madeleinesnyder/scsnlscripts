function lp = logData(P,D)
if isa(P,'th_GP'),
    lp = P.logp(D); 
    return;
end

if P.type == 1, 
    lp = limvnpdf(D, P.mu,P.S_inv);
elseif P.type == 2, 
    lp = limvnpdf(D, P.mu,P.S_inv);   
elseif P.type == 5, 
    lp = limvnpdf(D,P.mu,P.S_inv);
else
    lp = 0; 
    return;
    assert(false);
    lp = linormpdf(D,P.mu,P.S_inv);    
end

end