function pos=mlpseek(pid, offset, origin);
% pos=mlpseek(pid, offset [, origin])
%
% Set pipe position indicator. Only works for reads.
% Almost identical in behaviour to ftell, but backwards seeking is very
% inefficient (requires closing and reopening the pipe, then moving
% forward).

if (nargin<2)
    help(mfilename)
end

if (~exist('origin','var') | isempty(origin))
    origin=0;
end

[filename,permission,machineformat,pipecommand]=mlpopen(pid);
if (permission=='w')
    error('(mlpseek) Seeking only supported for reading pipes')
end

switch (origin)
  case {'cof',0}
    origin=0;
  case {'bof',-1}
    origin=-1;
  case {'eof',1}
    origin=1;
  otherwise
    error('(mlpseek) Incorrect offset parameter')
end

pos=mlpmex(pid,'seek',offset,origin);

% $$$ FSEEK Set file position indicator. 
% $$$     STATUS = FSEEK(FID, OFFSET, ORIGIN) repositions the file position
% $$$     indicator in the file associated with the given FID.  FSEEK sets the 
% $$$     position indicator to the byte with the specified OFFSET relative to 
% $$$     ORIGIN.
% $$$  
% $$$     FID is an integer file identifier obtained from FOPEN.
% $$$  
% $$$     OFFSET values are interpreted as follows:
% $$$         >= 0    Move position indicator OFFSET bytes after ORIGIN.
% $$$         < 0    Move position indicator OFFSET bytes before ORIGIN.
% $$$  
% $$$     ORIGIN values are interpreted as follows:
% $$$         'bof' or -1   Beginning of file
% $$$         'cof' or  0   Current position in file
% $$$         'eof' or  1   End of file
% $$$  
% $$$     STATUS is 0 on success and -1 on failure.  If an error occurs, use
% $$$     FERROR to get more information.
% $$$  
% $$$     Example:
% $$$  
% $$$         fseek(fid,0,-1)
% $$$  
% $$$     "rewinds" the file.
% $$$  
