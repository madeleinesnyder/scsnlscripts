function varargout=zseek(varargin)
% varargout=zseek(zid,varargin)
%
% Gateway function to seek through a normal or compressed files (pipes).
% Calls fseek or pseek, depending on zid; see these functions for 
% detailed documentation.
% 

if (nargin<1)
    help(mfilename)
    return
end

zid=varargin{1};
if (zid<0)
    if (nargout>0)
        varargout{:}=mlpseek(varargin{:});
    else
        mlpseek(varargin{:});
    end
else
    if (nargout>0)
        varargout{:}=fseek(varargin{:});
    else
        fseek(varargin{:});
    end
end
