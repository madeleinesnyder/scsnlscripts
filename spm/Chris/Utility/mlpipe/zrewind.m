function varargout=zrewind(varargin)
% varargout=zrewind(zid,varargin)
%
% Gateway function to rewind a normal or compressed files (pipes).
% Calls frewind or prewind, depending on zid; see these functions for 
% detailed documentation.
% 

if (nargin<1)
    help(mfilename)
    return
end

zid=varargin{1};
if (zid<0)
    if (nargout>0)
        varargout{:}=mlprewind(varargin{:});
    else
        mlprewind(varargin{:});
    end
else
    if (nargout>0)
        varargout{:}=frewind(varargin{:});
    else
        frewind(varargin{:});
    end
end
