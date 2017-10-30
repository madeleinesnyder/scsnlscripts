function outlist = readlist(inlist)

% Readout groups and subjects in inlist
%
% inlist can be the following formats:
% 1. {'01-01-01', '02-02-02',....}:
%   one group, multiple subjects
% 2. '01-01-01'                                
%   one group, single subject
% 3. 'list.txt'                                
%   one group, multiple subjects
% 4. {{'01-01-01',...}, {'02-02-02',...}, ...}
%   multiple groups
% 5. {'list1.txt', 'list2.txt', ...} 
%   multiple groups
%
% outlist{group_number}{subject_id}
%__________________________________________________________________________
% 2009-2010 The Stanford Cognitive and Systems Neuroscience Laboratory
% Tianwen Chen, 04/09/2010

scsnl_id = '$readlist.m v1$';

sprintf('%s \n', scsnl_id);

if iscell(inlist)
  if iscell(inlist{1})
    outlist = inlist;
  else
    [pathstr, name, ext] = fileparts(inlist{1});
    if strcmp(ext,'.txt')
      numgroup = length(inlist);
      outlist = cell(numgroup,1);
      for gcnt = 1:numgroup
        outlist{gcnt} = getitem(inlist{gcnt});
      end
    else
      outlist{1} = inlist;
    end
  end
else
  outlist = cell(1,1);
  [pathstr, name, ext] = fileparts(inlist);
  if strcmp(ext,'.txt')
    outlist{1} = getitem(inlist);
  else
    outlist{1} = {inlist};
  end
end

end

function slist = getitem(filename)

if ~exist(filename, 'file')
  error('Cannot find file: %s \n', filename);
else 
  fid = fopen(filename);
  cnt = 1;
  while ~feof(fid)
    fstr = fgetl(fid);
    str  = strtrim(fstr);
    if ~isempty(str)
      slist{cnt} = str;
      cnt = cnt + 1;
    end
  end
  fclose(fid);
end

end