function lp = logPrior(P)
if length(P) > 1,
    lp = 0;
    for k=1:length(P), lp = lp + logPrior(P(k)); end
    return;
end
if P(1).type == 1,
   lp = logNWpdf(P);
elseif P.type == 2, %GP.
    lp = lgampdf(P.l) + lgampdf(P.sigf) + lgampdf(P.sign);
elseif P.type == 3, % Gamma matrix.
    lp = lgampdf(P.S_inv,P.alpha,P.beta);
    lp = sum(lp(:));
elseif P.type == 5,
    %% do stuff here
    % student half t prior.
    lp = sum(lhalftpdf(P.sigma,P.nu,P.A)) + sum(limvnpdf(P.mu', P.p_mu.mu0, P.p_mu.S0_inv));
else
    assert(false);
end
end