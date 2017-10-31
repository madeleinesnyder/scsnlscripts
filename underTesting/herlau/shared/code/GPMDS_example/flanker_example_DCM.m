function flanker_example_DCM()
% See flanker_example_GPMDS.m for (a few) comments
close all; clc
cdir = fileparts(mfilename('fullpath'));
cd(cdir);
%WIN = false; 
WIN = true;
addpath(cdir,'../'); 
mdsd = fullfile(cdir,'../');  
paths = {mdsd, fullfile(mdsd,'/GPMDS'),fullfile(cdir,'../DCMwrap'),fullfile(cdir,'../MDS')  };
fout = sprintf('flanker_example_DCM'); if WIN, fout=[fout,'WIN']; end
%<grid_rm>
if ~WIN, 
    addpath(fullfile(mdsd, '../farm'));
    o = struct();
    o.ON_GRID = true; 
    o.ncores = 1;
    o.walltime = 10; 
    o.maxque = 64; 
    o.paths = paths; 
    
    tic();
    rs = autogrid(mfilename(),'rs',o);
    toc();
    save(fout,'rs');
    pause(4);
    return;
end
%</grid_rm>
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
% <grid_before_loop> 
for s = 1:SUBJECTS
    for m=1:length(methods)
        for ses=1:SESSIONS
% <grid_within_loop>
%% <grid_fun>  
[~,ps] = fMRI_data(DSET, 5, 3);
task = 1; 
p = ps{ses,task};
p.s = p.s(s); 
opts.random_init = false; 
p.opts = opts;
  
if strcmp(methods{m},'DCM'),
    [stats,PP] = DCM_matlab(p,opts);
elseif strcmp(methods{m},'MDS')
    [stats,PP] = MDS_matlab(p,opts);   
else
    assert(false);
end
rs{s,m,ses}.stats = stats;
% </grid_fun>
        
        end
    end
end
% <grid_after_loop>
end