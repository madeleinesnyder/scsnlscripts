% MLF2PDF prints a PDF figure for inclusion in LaTeX documents, with
%    LaTeX fonts. It creates a complete pdf-file with annotation of the 
%    figure (titles, labels and texts) and graphics (lines,arrows,markers,
%    ...). This function uses the matlabfrag function and makes simply a
%    concatenation of the text part and the graphical part in a single PDF
%    file. It requires:
%
%                  MATLAB        |      OTHER
%               -----------------|----------------------------
%                 matlabfrag.m   |      pdflatex
%
%    
%    USAGE:
%    ------------------------
%
%    mlf2pdf creates the PDF file of the current graphical figure (CGF) 
%    named LaTeXfile.pdf
% 
%    mlf2pdf(H) creates the PDF file from the graphical figure with 
%    handle H, named LaTeXfile.pdf
%    
%    mlf2pdf(H,FILENAME) creates the PDF file from the graphical figure 
%    with handle H, named FILENAME.PDF (FILENAME is a character array of
%    the filename, with or without the extension .PDF)
%   
%    mlf2pdf(H,FILENAME,PACKAGES) to adding
%    packages such as fonts, mathfont,...
%    ex: PACKAGES = 'amssymb, times'
% 
function mlf2pdf(h,FileName,Packages,opts)
if nargin < 4, opts = struct(); end

if ~verLessThan('matlab', '9.0.1')
    disp('new matlab version, figure will look bad')
    version()
    set(h,'Units','inches');
    screenposition = get(gcf,'Position');
    set(h,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    print(FileName, '-dpdf','-painters');
    pause(0.1);
    docrop(FileName);
    
    return;
end
dopts.PAD_x = 0;
dopts.PAD_y = 0;
dopts.crop = true;
o2 = opts;
dopts.LooseInsetScale = 1.8;
set(h,'Units','Centimeters');
v = get(h,'Position');
if v(3) < 10 && v(3) > 6,
    dopts.PositionFudge = 1.1;
elseif v(3) < 6,
    dopts.PositionFudge = 1.1;    
else
    dopts.PositionFudge = 1;
end
opts = setstructfields(dopts,opts);

xl = get(gca,'XLabel');
has_xl = length(xl.String) > 0;

yl = get(gca,'YLabel');
has_yl = length(yl.String) > 0;

% set(h,'Units','Centimeters','Position',v .* [1 1 dopts.PositionFudge dopts.PositionFudge])

% padx = [0 0];
% pady = [0 0];

% how about this: Assume the axis scale is okay, then pad by altering the
% outer/inner locations?
%%
PAD = 1; % padding in centimeters.

set(gca,'Units','Centimeters');
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 

axTI = ax.TightInset;
axPos = ax.Position;
axOPos = ax.OuterPosition;

fPos = get(h,'Position');
% [fPos ; axTI ; axPos ; axOPos]
 
hl = legend;
34;  
set(hl,'Units','Centimeters');
lPos = get(hl,'Position');
%  45     
%%
set(h,'Position',fPos + [0 0 2*opts.PAD_x 2*opts.PAD_y]); %DONE
% hchld=get(gcf,'children'); 
% for j=1:length(hchld),
%     if length(h)>1, % subplots,    
%     ax = gca;    
ax.Position = ax.Position + 0*[opts.PAD_x opts.PAD_y 0 0];
if ~isempty(hl) && false,
    
    lPos = lPos + [opts.PAD_x opts.PAD_y 0 0];
    set(hl,'Position',lPos)
end 

%     end
% end
% 
% set(gca(h), 'LooseInset', v*opts.LooseInsetScale);

%%
if nargin < 2
    FileName='LaTeXfile';
end

if nargin < 1
    h=gcf;
end

TempName = strcat('TEMP',num2str(round(rand*10000))); %Generate random file name

if ~exist('matlabfrag.m')
    disp('MatLabFrag M-file does not exist or is not in the MATLAB''s search path.');
    disp('This file can be downloaded at: http://www.mathworks.com/matlabcentral/fileexchange');
    disp('                                                         Try again...');
    return;
end

matlabfrag(TempName,'handle',h);  %call matlabfrag to export figure to .eps and .tex file.

%-------------------------------------------------------
% Temporary LaTeX file
%-------------------------------------------------------
fid = fopen(strcat(TempName,'2.tex'),'w');

fprintf(fid,'\\documentclass[11pt, oneside]{article}\n');
fprintf(fid,'\\usepackage{graphicx}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\usepackage[T1]{fontenc}\n');
fprintf(fid,'\\usepackage[latin1]{inputenc}\n');
fprintf(fid,'\\usepackage{ae}\n');
fprintf(fid,'\\usepackage{psfrag}\n');
fprintf(fid,'\\usepackage{color}\n');
fprintf(fid,'\\usepackage{pstool}\n');

if nargin>=3 && ~isempty(Packages),
    fprintf(fid,'\\usepackage{%s}\n',Packages); % Suplementary packages
end
fprintf(fid,'\\pagestyle{empty}\n');
fprintf(fid,' \n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'    \\begin{figure}\n');
fprintf(fid,'        \\centering\n');
fprintf(fid,'        \\psfragfig{%s}\n',TempName);
fprintf(fid,'    \\end{figure}\n');
fprintf(fid,' \n');
fprintf(fid,'\\end{document}\n');
fclose(fid);

%-------------------------------------------------------
% LaTeX Command
%-------------------------------------------------------
Str=sprintf('pdflatex -shell-escape --src -interaction=nonstopmode %s2.tex',TempName);
disp(sprintf('\n[LaTeX Command] %s',Str));
[hdos,wdos]=system(Str);
pause(0.1);
if hdos ~=0
   if isunix==0
     dos(sprintf('del %s*',TempName));  
   else
     unix(sprintf('rm %s*',TempName));
   end
   error('Error %d -- LATEX:\n%s',hdos ,wdos);
   return;
end


%-------------------------------------------------------
% Rename TempFile to FileName and delete FileName before (if it exists)
%-------------------------------------------------------
if isunix==0
%    delete(FileName);   
   movefile([TempName,'.pdf'],[FileName,'.pdf'])
%    dos(sprintf('del ''%s.pdf''',FileName));
%    dos(sprintf('ren %s.pdf %s.pdf',TempName,)); 
else 
   unix(sprintf('rm %s.pdf',FileName));
   unix(sprintf('mv %s.pdf %s.pdf',TempName,FileName));
end
if opts.crop,
    docrop(FileName);

end
%system(sprintf('pdftops -eps %s.pdf',FileName))

%-------------------------------------------------------
% Success
%------------------------------------------------------- 
  disp(sprintf('... OK!\nPDF file [%s.pdf] has been created in the current directory\n',FileName));
 %-------------------------------------------------------
% Delete all the temporary files
%-------------------------------------------------------
if isunix==0
  dos(sprintf('del %s*',TempName));
  dos(sprintf('del TEMP*.tex'));  
else
  unix(sprintf('rm %s*',TempName));
end
return;
end
function docrop(FileName)
    s = [FileName,'.pdf'];
    
    j = system(sprintf('pdfcrop %s %s',s,s));
    if j == 127
        %%
        %clc
        PATH = getenv('PATH');
        setenv('PATH', [PATH ':/Library/TeX/texbin:/usr/local/bin']);
        %unix('setenv PATH "/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/usr/local/bin"') 
        j = system(sprintf('pdfcrop %s %s',s,s));
        % /Library/TeX/texbin/pdfcrop
    end
end