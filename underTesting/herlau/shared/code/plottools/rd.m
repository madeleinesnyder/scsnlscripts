function M = rd(M)
if iscell(M),
    for j=1:numel(M),
        M{j} = rd(M{j});
    end
    return;
end
for j=1:size(M,3),
    M(:,:,j) = M(:,:,j)-diag(diag(M(:,:,j)));
end
end