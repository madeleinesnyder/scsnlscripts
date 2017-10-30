function TR = GetTR

files = dir('E*');

% If no E file exists, assume a TR of 2.
if(length(files) == 0)
  TR = 2;
  return;
end

fid   = fopen(files(1).name);

while 1
    tline = fgetl(fid);
    if ~ischar(tline),   break,   end
    [a b c d e] = regexpi(tline, '^TR[\s\t]*=[\s\t]*(\d+)');
    if(~isempty(a))
      TR = str2double(e{1})/1000;
      break;
    end
    
end

if fid > 0; fclose(fid); end

end


