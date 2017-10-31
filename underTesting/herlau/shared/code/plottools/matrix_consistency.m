function [val,M,valJ] = matrix_consistency(C)
C = C(:);
[n,~,J] = size(C{1});
%S = length(C);
%n = size(C{1},1);
II = find( repmat(1-eye(n), [1, 1, J]));
%%
for i=1:length(C),
   for j=1:length(C),
       M(i,j) = corr(C{i}(II), C{j}(II), 'type','Kendall');
   end
end
      %%
I = triu(ones(length(M)),1);
val = mean(M(I>0));

if J == 1, valJ = val; return; end
for j=1:J,
    C2 = cell(size(C));
    for i=1:length(C),
       C2{i} = C{i}(:,:,j);
       
    end
    valJ(j) = matrix_consistency(C2);
end

end