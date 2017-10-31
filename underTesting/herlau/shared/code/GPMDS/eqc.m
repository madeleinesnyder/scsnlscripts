function b = eqc(X,Y) % check if X approximately equals Y with some error checking. 
X = X(:); Y = Y(:);
I = ~isfinite(X(:));
b = true; 
if any(I), b = all( X(I) == Y(I) ); return; end
if any(~I),
    dm = abs(X(~I)-Y(~I));
    dp = 1 + abs(X(~I))+abs(Y(~I));
    b = max(dm ./ dp)  < 10^-8; 
    return;
end
end