function varargout=zopen(varargin)
% varargout=zopen(varargin)
%
% Gateway function to open normal or compressed files (pipes).
% Checks the filename extension - if it is one of .gz, .Z or .zip, 
% calls mlpopen(), otherwise calls fopen(); see these functions for 
% detailed documentation.
% 
% NOTE: zopen/zclose and related functions use a special file identifier,
% called a zid, which is positive for normal files and negative for pipes.
% Use the zid for any calls to zclose, zwrite, zread etc.
% For normal files, the zid is identical to the fid returned by fopen (+ve).
% For compressed files, the zid is identical to the pid returned by mlpopen (-ve).
% 
% zids=zopen('all') returns a vector of the zids of all open files and pipes.
%
% EXAMPLE
%  z1=zopen('mycompressedfile.gz','w')
%  z2=zopen('mynormalfile','w')
%  c=zwrite(z1,rand(100),'float')
%  pos=ztell(z2)
%  c=zwrite(z2,0:100,'uint8')
%  pos=ztell(z2)
%  zids=zopen('all')
%  zclose(z1)
%  zclose(z2)
%  z1=zopen('mycompressedfile.gz','r')
%  z2=zopen('mynormalfile','r')
%  zseek(z2,50,0)
%  [d1,c1]=zread(z1,0,'float');
%  [d2,c2]=zread(z2);
%  zclose('all')

if (nargout<1 | nargin<1)
    help(mfilename)
    return
end

varargout=cell(nargout,1);

if (nargin==1)
    % return info about pipe or file
    if (isnumeric(varargin{1}))
        if (varargin{1}<0)
            if (nargout<4)
                help('mlpopen')
                return
            end
            [varargout{1},varargout{2},varargout{3},varargout{4}]=mlpopen(varargin{1});
        else
            if (nargout<4)
                help('fopen')
                return
            end
            [varargout{1},varargout{2},varargout{3},varargout{4}]=fopen(varargin{1});        
        end
    elseif (ischar(varargin{1}) & strcmp(varargin{1},'all'))
        p=mlpopen('all');
        f=fopen('all');
        varargout{1}=[f p];
    end
    return
else
    fname=varargin{1};
    if (nargin>1 && strcmp(varargin{2},'r') && ~exist(fname,'file'))
        error('(zopen) File not found')
    end
    %[pathstr,bname,ext,ver]=fileparts(fname);
    [pathstr,bname,ext]=fileparts(fname);
    ver = -1; %JK
    if (strcmp(ext,'.gz') | strcmp(ext,'.Z') | strcmp(ext,'.zip'))
        varargout{:}=mlpopen(varargin{:});
    else
        varargout{:}=fopen(varargin{:});
    end
end

