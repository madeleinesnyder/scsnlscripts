function flanker_example_DCM()
% See flanker_example_GPMDS.m for (a few) comments
close all; clc
cdir = fileparts(mfilename('fullpath'));
cd(cdir);
WIN = false; 
%WIN = true;
addpath(cdir,'../'); 
mdsd = fullfile(cdir,'../');  
paths = {mdsd, fullfile(mdsd,'/GPMDS'),fullfile(cdir,'../DCMwrap'),fullfile(cdir,'../MDS')  };
fout = sprintf('flanker_example_DCM'); if WIN, fout=[fout,'WIN']; end
cd(cdir)
addpath(paths{:});
DSET = 'flanker';
opts.TT = 200; 
%%
[p,ps] = fMRI_data(DSET); SUBJECTS = length(p.s);
[SESSIONS,TASKS] = size(ps);
p = []; ps = [];
%%
save([fout '_stuff']);
methods = {'DCM', 'MDS'};
rs = cell(SUBJECTS,length(methods),SESSIONS);
pari = 1; pars = cell(size(rs));
S = functions(str2func(mfilename()));
[f,d] = fileparts(S.file);
addpath(f);

for s = 1:SUBJECTS
    for m=1:length(methods)
        for ses=1:SESSIONS
x = whos();
pars{pari} = struct(); 
for GLi=1:length(x),
    if any(strcmp(x(GLi).name, {'x','pars','ans'})), continue; end
    pars{pari}.(x(GLi).name) = eval(x(GLi).name);
end
pari = pari+1;

        
        end
    end
end

 grid_farm_prepare('flanker_example_DCM',pars);

end
