function delete_dir_inside (dir_list)

current_dir = pwd;

dlist = ReadList(dir_list);

numfile = length(dlist);

for i = 1:numfile
  cd(dlist{i});
  unix('rm -rf *');
end

cd(current_dir);

end