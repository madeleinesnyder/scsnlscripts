function y = lngengamma(M,alpha)
y = log(pi) * (M*(M-1)/4) + sum(gammaln( alpha + (1 - [1:M] )/2 )  );
end