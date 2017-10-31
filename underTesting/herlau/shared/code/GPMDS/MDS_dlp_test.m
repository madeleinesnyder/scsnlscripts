function MDS_dlp_test(pn,po,xn,xo,mu,S_inv)
dlp0 = MDS_logp(pn) - MDS_logp(po);      
dlp1 = -1/2 * (xn - mu)' * S_inv * (xn - mu)  + 1/2 * (xo - mu)' * S_inv * (xo - mu);

dd = abs(dlp0-dlp1)/abs(dlp0+dlp1);

if dd > 10^-4,
    disp([dlp0,dlp1]);
    disp(dd)     
    assert(false);
end  
% disp('passed');   


end