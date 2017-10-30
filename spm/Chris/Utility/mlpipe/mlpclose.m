function mlpclose(pid)
% mlpclose(pid)
% mlpclose('all')

if (nargin==0)
    help(mfilename)
    return
end

if (ischar(pid))
    if (strcmp(pid,'all'))
        mlpmex(0,'close');
    else
        help(mfilename)
        return
    end
else
    mlpmex(pid,'close');
end
