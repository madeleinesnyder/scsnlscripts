function matlab_jinjafy(fin,fout)
if nargin < 1, KendallMatrix(); return; end
CDIR = pwd();
jindir = fileparts(mfilename('fullpath'));

%% jinjafy the text file. Recall you can use load commands within the file. 
[fin_dir, fin_name,ex] = fileparts(fin);
if nargout < 2, 
    if any(fin_name == '_'),
        fout_name = [fin(1:find(fin_name=='_',1,'last')-1), ex];
    else
        assert(false);
    end    
    fin_dir = fullfile(CDIR, fin_dir); 
    fout_dir = fin_dir;
else
    [fout_dir, fout_name, ex] = fileparts(fout);
    fout_name = [fout_name, ex];
    fout_dir = fullfile(CDIR, fout_dir); 
end
jinja_py = fullfile(jindir, [fin_name, '_jinja.py']);
fin_name = [fin_name, ex];
%%
disp('Jinjafying text file...');
disp(['READING: ' fullfile(fin_dir,fin_name)]);
disp(['WRITING: ' fullfile(fout_dir,fout_name)]);

%%
%jindir
cd(jindir);
f = fopen(jinja_py,'w');
fprintf(f, 'from jinja_weaver import jinja_process\n');
fprintf(f, 'fin_dir = "%s"\n',fin_dir);
fprintf(f, 'fin_name = "%s"\n',fin_name);
fprintf(f, 'fout_dir = "%s"\n',fout_dir);
fprintf(f, 'fout_name = "%s"\n',fout_name);
fprintf(f, 'jinja_process(fin_dir,fin_name,fout_dir,fout_name)\n');
%%

fclose(f);
pause(0.1); 
[~,fn,fex] = fileparts(jinja_py);
system(['python ' fn fex]);
cd(CDIR); 
end