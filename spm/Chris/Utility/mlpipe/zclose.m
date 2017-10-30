function varargout=zclose(varargin)
% varargout=zclose(zid)
%
% Gateway function to close normal or compressed files (pipes).
% 
% zclose(pid)                  close pipe or file
% zclose([pid1 pid2 zid1 ..])  close multiple files and/or pipes
% zclose('all')                close all pipes and files

if (nargin<1)
    help(mfilename)
    return
end

if (ischar(varargin{1}))
    if (strcmp(varargin{1},'all'))
        fclose('all');
        mlpclose('all');
    else
        help(mfilename)
        error('(zclose) Unknown syntax')
    end
elseif (isnumeric(varargin{1}))
    for n=1:length(varargin{1})
        zid=varargin{1}(n);
        if (zid<0)
            mlpclose(zid);
        else
            fclose(zid);
        end
    end
else
    help(mfilename)
    error('(zclose) Unknown syntax')
end
