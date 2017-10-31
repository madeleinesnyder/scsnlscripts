function par = flanker_example_DCM_grid_fun(par,pig)

DSET = par.DSET;             
S = par.S;                   
SESSIONS = par.SESSIONS;     
SUBJECTS = par.SUBJECTS;     
TASKS = par.TASKS;           
WIN = par.WIN;               
cdir = par.cdir;             
d = par.d;                   
f = par.f;                   
fout = par.fout;             
m = par.m;                   
mdsd = par.mdsd;             
methods = par.methods;       
opts = par.opts;             
p = par.p;                   
pari = par.pari;             
paths = par.paths;           
ps = par.ps;                 
rs = par.rs;                 
s = par.s;                   
ses = par.ses;               
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
par.rs = rs;

 end 
