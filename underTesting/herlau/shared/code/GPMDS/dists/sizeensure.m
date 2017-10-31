for j=1:length(Svars),
    SS{j} = size(eval(Svars{j}));        
end
l = max(cellfun(@length, SS));
M = zeros(length(SS), l);
for j=1:length(Svars),
    M(j,:) = SS{j};    
end
MM = repmat(max(M, [],1), length(SS), 1);
for j=1:length(Svars),
    dm = repmat(eval(Svars{j}), MM(j,:) ./ M(j,:)); 
   eval([Svars{j} , '= dm;']);    
end