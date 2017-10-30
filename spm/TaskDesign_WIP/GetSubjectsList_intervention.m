function [scanNumb, subNumb] = GetSubjectsList(filename)

% fid = fopen(filename);
% str = fgetl(fid);
% cnt = 1;
% 
% while(str ~= -1)
%   subjects{cnt} = str;
%   str = fgetl(fid);
%   cnt = cnt + 1;
% end
% 
% fclose(fid);

fid = fopen(filename);
str = fgetl(fid);
cnt = 1;

ica = [];

while(str ~= -1)
  inds      = regexpi(str, '[\t,]');
  if(~isempty(inds))
    scanNumb{cnt} = str(1:inds(1)-1);
    inds = [inds, length(str) + 1];
    subNumb(:, cnt) = str2num(str(inds(1)+1:end))';
  else
    scanNum{cnt} = str;    
  end
  str = fgetl(fid);
  cnt = cnt + 1;
end


fclose(fid);