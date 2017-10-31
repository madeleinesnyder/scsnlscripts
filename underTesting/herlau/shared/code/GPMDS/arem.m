function s2 = arem(s,fld)
s2 = cell(size(s)); 
for i=1:numel(s),
    s2{i} = s(i);
    if isfield(s(i), fld), 
        s2{i}= rmfield(s2{i}, fld);
    end
end
s2 = reshape([s2{:}], size(s));

end