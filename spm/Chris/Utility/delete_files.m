function delete_files (file_list)

flist = ReadList(file_list);

numfile = length(flist);

for i = 1:numfile
  unix(sprintf('rm -rf %s', flist{i}));
end

end