function [H,W] = getHW(n)
W = ceil(sqrt(n));
H = ceil(n/W);
if nargout == 1, 
    H = [H,W];
end
end