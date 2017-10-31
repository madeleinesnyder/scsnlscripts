function GPMDS_Example1()
% This example evaluate the sampler on a single rat
cdir = fileparts(mfilename('fullpath'));
cd(cdir);
addpath(fullfile(cdir,'../'));
addpath(fullfile(cdir,'../MDS'));
addpath(fullfile(cdir,'../DCMwrap'));

p = fMRI_data('rats');
p.s = p.s(1); % only use first rat (alternative: p.s = p.s(2:3) to use rat 2&3). 
opts.T = 200; % run for 1000 iterations. 
opts.plot = true; 
opts.stat_iter = 10; 
%[stats,p2] = GPMDS_matlab(p,opts);
[a,b] = DCM_matlab(p);
[c,d] = MDS_matlab(p);

end