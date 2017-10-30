function move_files (from_dir, to_dir, list)

file_list = ReadList(list);
numfile = length(file_list);

for i = 1:numfile
  org_dir = fullfile(from_dir, file_list{i});
  dest_dir = to_dir;
  unix(sprintf('cp -a %s %s', org_dir, dest_dir));
  dest_dir = fullfile(to_dir, file_list{i});
  unix(sprintf('diff -rf %s %s > %s', org_dir, dest_dir, ...
    ['dir_', num2str(i),'.txt']));
end

end