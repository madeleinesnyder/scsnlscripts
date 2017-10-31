function y = lhalfcachypdf(x,mu,gamma)
y = lhalftpdf(x - mu, 1, gamma);
end
