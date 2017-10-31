function flanker_example_GPMDS()
close all; clc
cdir = fileparts(mfilename('fullpath'));
cd(cdir);
%WIN = true; % set to true to run on local machine. 
WIN = false; 
addpath(cdir,'../'); 
mdsd = fullfile(cdir,'../');  
paths = {mdsd, fullfile(mdsd,'/GPMDS') };
% outputfile to save to. 
fout = [sprintf('flanker_example_GPMDS')]; if WIN, fout=[fout,'WIN']; end
%<grid_rm> % These tags are there to parse file. This block has to be
%removed to prevent infinite recursion when run on server. 
if ~WIN,
    % This is the code that execute on server
    addpath(fullfile(mdsd, '../farm'));
    o = struct();
    o.ON_GRID = true; 
    o.ncores = 1; % cores PER JOB. 
    o.walltime = 5; 
    o.maxque = 64;  % maximum qued jobs. 
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
TT = 200; % number of iterations. Start with about 2000 or more.  
opts.stat_iter = round(TT/200); 
if WIN
    TT = 1000; 
    opts.stat_iter = 5; 
    opts.plot =true; 
end

[p,ps] = fMRI_data(DSET); SUBJECTS = length(p.s); p = [];
[SESSIONS,TASKS] = size(ps);
p = []; ps = []; % remember to remove variables to prevent too much storage. 
CONT = false; % Continue MCMC simulation from last stop. Must initially be false. 
save([fout '_stuff']);
% must be pre-allocated. 
rs = cell(SUBJECTS,SESSIONS);     
% <grid_before_loop> 
for s = 1:SUBJECTS
    for ses=1:SESSIONS            
% <grid_within_loop>
%% <grid_fun>  
[~,ps] = fMRI_data(DSET, 5, 3);
task = 1; 
p = ps{ses,task};
p.s = p.s(s); 
opts.T = TT;
II = [s,ses]; % identifier for stopcon wrap. 
% Stopcon wrap can be removed. It is a special function that allows sampler
% to continue from last stop by naming resource files appropriately. 
stats = GPMDS_stopcon_wrap(p,opts,fout,II,@(p,opts)GPMDS_matlab(p,opts),CONT);
rs{s, ses}.stats = stats; % what is actually to be saved. 
% </grid_fun>
    end
end    
% <grid_after_loop>
end