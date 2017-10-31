function s = MDS_Prnd(P,M,T) % random element of whatever that structure represents.
if P.type == 1, % NW.    
    s = imvnrnd(P.mu', P.S_inv,T)';     
elseif P.type == 2, 
    %% 
    s = imvnrnd(P.mu', P.S_inv,M);
else    
    assert(false);
end


end