function file_deletion

flist = ReadList('/home/tianwenc/file_for_deletion.txt');
pfolder = '/mnt/plum3_share3/integers_data';

nfile = length(flist);

for i = 1:nfile
  file_name = fullfile(pfolder, flist{i});
  unix(sprintf('/bin/rm -rf %s', file_name));
end

end