function varargout=zread(varargin)
% varargout=zread(zid,varargin)
%
% Gateway function to read from normal or compressed files (pipes).
% Calls fread or pread, depending on zid; see these functions for 
% detailed documentation.
% 

if (nargin<1 | nargout<1)
    help(mfilename)
    return
end

zid=varargin{1};
if (nargout<2)
    if (zid<0)
        varargout{1}=mlpread(varargin{:});
    else
        varargout{1}=fread(varargin{:});
    end
else
    if (zid<0)
        [varargout{1},varargout{2}]=mlpread(varargin{:});
    else
        [varargout{1},varargout{2}]=fread(varargin{:});
    end    
end