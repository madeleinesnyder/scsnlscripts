function y = halfcachyrnd(mu,gamma)
%%
ov = 1; 
Svars = {'ov', 'mu','gamma'}; 
sizeensure
%%
y = halftrnd(ov,gamma) + mu;
end