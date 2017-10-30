function varargout=ztell(varargin)
% varargout=ztell(zid,varargin)
%
% Gateway function to det file position for normal or compressed files (pipes).
% Calls ftell or ptell, depending on zid; see these functions for 
% detailed documentation.
% 

if (nargin<1 | nargout<1)
    help(mfilename)
    return
end

zid=varargin{1};
if (zid<0)
    varargout{:}=mlptell(varargin{:});
else
    varargout{:}=ftell(varargin{:});
end
