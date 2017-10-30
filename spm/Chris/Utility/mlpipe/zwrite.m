function varargout=zwrite(varargin)
% varargout=zwrite(zid,varargin)
%
% Gateway function to write to normal or compressed files (pipes).
% Calls fwrite or pwrite, depending on zid; see these functions for 
% detailed documentation.
% 

if (nargin<1)
    help(mfilename)
    return
end

zid=varargin{1};
if (zid<0)
    if (nargout>0)
        varargout{:}=mlpwrite(varargin{:});
    else
        mlpwrite(varargin{:});
    end
else
    if (nargout>0)
        varargout{:}=fwrite(varargin{:});
    else
        fwrite(varargin{:});
    end
end
