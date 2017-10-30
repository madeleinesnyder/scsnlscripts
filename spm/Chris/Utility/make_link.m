function make_link (org_path, dest_path, list)

file_list = ReadList(list);
numfile = length(file_list);
current_dir = pwd;
cd(org_path);
for i = 1:numfile
  org_dir = fullfile(org_path, file_list{i});
  dest_dir = fullfile(dest_path, file_list{i});
  unix(sprintf('/bin/rm -rf %s', org_dir));
  unix(sprintf('ln -sT %s %s', dest_dir, file_list{i}));
end

cd(current_dir);

end