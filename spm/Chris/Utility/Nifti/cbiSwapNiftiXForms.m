function h=cbiSwapNiftiXForms(filename,newfilename)
% h=cbiSwapNiftiXForms(filename [,newfilename])
% 
% Swaps the qform and sform info of a file and saves the result.
% If a newfilename is specified, the result is saved to this file.
% 

if (nargin<1)
    help(mfilename)
    return
end

[d,h]=cbiReadNifti(filename);


qf=h.qform44
qfc=h.qform_code;
sf=h.sform44
sfc=h.sform_code;

h.qform44=sf;
h.qform_code=sfc;
h.sform44=qf;
h.sform_code=qfc;

if (nargin==2)
    filename=newfilename;
end

cbiWriteNifti(filename,d,h)


