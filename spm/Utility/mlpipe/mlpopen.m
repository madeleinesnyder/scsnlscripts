function varargout=mlpopen(varargin)
% pid=mlpopen(filename,permission [,machineformat,pipecommand] )
% [filename,permission,machineformat,pipecommand]=mlpopen(pid)
% pids=mlpopen('all')

if (nargout<1)
    help(mfilename)
    return
end

if (nargin==1)
    if (ischar(varargin{1}) & strcmp(varargin{1},'all'))
        varargout{1}=mlpmex(0,'open');
        return
        
    elseif (isnumeric(varargin{1}))
        % [filename,permission,machineformat,pipecommand]=mlpopen(pid)
        [filename,permission,byteswap,pipecommand]=mlpmex(varargin{1}(1),'open');
        if (nargout>0)
            varargout{1}=filename;
        end
        if (nargout>1)
            varargout{2}=permission;
        end
        if (nargout>2)
            % convert byteswap to endian
            [c,maxsize,endian]=computer;
            if (byteswap) 
                if (lower(endian)=='l')
                    varargout{3}='b';
                else
                    varargout{3}='l';
                end
            else
                varargout{3}=lower(endian);
            end
        end
        if (nargout>3)
            varargout{4}=pipecommand;
        end
        return
        
    end
elseif (nargin>1)
    % pid=mlpopen(filename,permission [,machineformat,pipecommand] )
    filename=varargin{1};
    permission=varargin{2};
    if (nargin>2)
        machineformat=varargin{3};
    else
        machineformat='native';
    end
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
        
    if (nargin>3)
        pipecommand=varargin{4};
    else
        [pstr,name,ext,ver]=fileparts(filename);
        if (permission=='r')
            switch (ext)
              case '.Z'
                pipecommand='uncompress -c ';
              case '.gz'
                pipecommand='gunzip -c ';
              case '.zip'
                pipecommand='unzip -c ';
              otherwise
                error('popen(): Unknown file type - cannot open pipe!')
            end
        elseif (permission=='w')
            switch (ext)
              case '.Z'
                pipecommand='compress > ';
              case '.gz'
                pipecommand='gzip > ';
              case '.zip'
                pipecommand='zip > ';
              otherwise
                error('(mlpopen) Unknown file type - cannot open pipe!')
            end
        else
            error('(mlpopen) permission must be one of ''r'' or ''w''');
        end            
    end    
    pid=mlpmex(0,'open',filename,permission,byteswap,pipecommand);
    if (nargout>0)
        varargout{1}=pid;
    else
        disp(pid)
    end
    return
    
else
    help(mfilename)
    return
end
    