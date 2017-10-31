function kendall_test()
M = 5; 
n = 16 ; 
T = 100;
for t=1:T,
    k(t,1) = dkendall(n,M);
end
%%
close all;
mean(k), std(k)

plot(k,'.');
boxplot(k);


%%
end
function v2 = dkendall(n,M)
K = M * M - M;
%%
sigs = zeros(n,K); 
J = ~eye(M);
for i=1:n,
    x = rand(M,M); 
    g = randi(M);
    x(:,g) = x(:,g)+1;
    
    [~,I] = sort(x(J(:)));
   sigs(i,:) = I;
end
for i=1:n,
    for j=1:n,
%         X(i,j) = corr(sigs(i,:)', sigs(j,:)', 'type', 'Kendall');
    end
end
% v = X(triu(ones(n),1)==1);
% v = mean(v(:));

X2 = corr(sigs', sigs', 'type', 'Kendall');
v2 = X2(triu(ones(n),1)==1);
v2 = mean(v2(:));

end