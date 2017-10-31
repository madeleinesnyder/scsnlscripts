function y = dirpdf(x,alpha)
% if x(end) == 0 && alpha(end) == 0, 
%     x = x(1:end-1);
%     alpha = alpha(1:end-1);
% end

y = gammaln(sum(alpha)) - sum(gammaln(alpha)) + sum(log(x) .* (alpha-1));

end