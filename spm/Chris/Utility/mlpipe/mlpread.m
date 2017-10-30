function varargout=mlpread(pid,size,precision,machineformat)
% [data [,count] ]=mlpread(pid [,size,precision,machineformat] )
 
if (nargout<1 | nargin<1)
    help(mfilename)
    return
end

if (~exist('size','var') || isempty(size))
    size=Inf;
end
if (~isfinite(size)) 
    size=0;
end

if (~exist('precision','var') || isempty(precision))
    precision='uint8';
end
if (~exist('machineformat','var') || isempty(machineformat))
    machineformat='native';
end

% convert precision to nbytes
switch (precision)
  case {'uint8','uchar'}
    datatype='uchar';
    nbytes=1;
  case {'int8','char'}
    datatype='char';
    nbytes=1;
  case {'uint16','ushort'}
    datatype='ushort';
    nbytes=2;
  case {'int16','short'}
    datatype='short';
    nbytes=2;
  case {'int32','int'}
    datatype='int';
    nbytes=4;
  case {'float32','float'}
    datatype='float';
    nbytes=4;
  case {'float64','double'}
    datatype='double';
    nbytes=8;    
  otherwise
    error('Unknown or unsupported precision.');
end

% determine byteswapping
[c,maxsize,endian]=computer;
switch (machineformat)
  case {'native','n'}
    byteswap=0;
  case {'ieee-le','l'}
    if (lower(endian)=='l')
        byteswap=0;
    else
        byteswap=1;
    end
  case {'ieee-be','b'}
    if (lower(endian)=='l')
        byteswap=1;
    else
        byteswap=0;
    end
  otherwise
    error('(mlpopen) Unsupported machine format!')
end

% call mlpmex
[data,count]=mlpmex(pid,'read',size,nbytes,datatype,byteswap);

% output

varargout{1}=data;
if (nargout==2)
    varargout{2}=count;
end
