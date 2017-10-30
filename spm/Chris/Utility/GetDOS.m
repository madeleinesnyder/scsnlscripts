function yearid = GetDOS

files = dir('E*');

% If no E file exists, assume a TR of 2.
if isempty(files)
  disp('No E file found');
  yearid = -1;
  return;
end

fid   = fopen(files(1).name);

while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    [a b c d e] = regexpi(tline, ...
                  '^date of scan[\s\t]*=[\s\t]*\d*/\d*/1(\d*)');
    if(~isempty(a))
      yearid = str2double(e{1});
      break;
    end
    
end

if fid > 0; fclose(fid); end

end