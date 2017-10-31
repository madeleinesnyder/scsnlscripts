function [p,psall] = fMRI_data(dset,varargin)
cdir = fileparts(mfilename('fullpath'));
if nargin < 1, dset = 'rats'; end 
psall = [];
addpath(fullfile(cdir,'../spm12'));
addpath(fullfile(cdir,'rats'));
addpathif(fullfile(cdir,'7_node_subjects/MDSScripts'));
addpathif(fullfile(cdir,'MDSOpenfMRIData'));
addpathif(fullfile(cdir,'dcm_bms'));
addpathif(fullfile(cdir,'attention'));
addpathif(fullfile(cdir,'SPMrest/spDCM'));
addpathif(fullfile(cdir,'MDSadata'));
addpathif(fullfile(cdir,'mmc'));
addpathif(fullfile(cdir,'Spike'));
addpathif(fullfile(cdir,'GPMDS/dists'));
addpathif(fullfile(cdir,'GPMDS/'));
addpathif(fullfile(cdir,'Ox_sims/data'));
addpathif(fullfile(cdir,'Luna'));
addpathif(fullfile(cdir,'Yale'));
addpathif(fullfile(cdir,'StanfordStopGo'));
addpathif(fullfile(cdir,'oddball'));
addpathif(fullfile(cdir,'hcp/MDSScripts'));


addpathif(fullfile(cdir,'flanker'));
addpathif(fullfile(cdir,'SSTxui'));
 
p = struct(); 
if strcmp(dset,'rats'),
    p = ratdat(varargin{:});
    p.s = p.s(1:3);
end
if strcmp(dset,'hcp'),
    p = HCPWM_data(varargin{:});
end
if strcmp(dset,'hcp2ses'),
    p = HCPWM_data(varargin{:});
    %%
    S = length(p.s);
    psall = cell(2,1); 
    I = 1:S <= S/2;
    psall{1} = p; 
    psall{1}.s = psall{1}.s(I); 
    psall{2} = p; 
    psall{2}.s = psall{2}.s(~I); 
    p = psall{1};
    
end

if strcmp(dset,'stopgo'),
    %% load stanford stopgo data.  
    if ~isempty(strfind(pwd(), '/Users'))
        dr = fullfile(cdir,'../../../Documents/hippa/StopGoStanford/');
    else
        dr = fullfile(cdir,'../../hippa/StopGoStanford/');
    end
    
    addpath(dr);
    p = stopgo_data(varargin{:});
end
if strcmp(dset,'open'),
    p = OpenfMRI_data();    
end 
if strcmp(dset,'bms'),
    p = bms_data(varargin{:});
end 
if strcmp(dset,'attention'),
    p = attention_data(varargin{:});
end 
if strcmp(dset,'spDCM'),
    p = spDCM_data(varargin{:});
end 
if strcmp(dset,'mmc1'),
    p = mmc_data(1);
end 
if strcmp(dset,'mmc2'),
    p = mmc_data(2);
end 
if strcmp(dset,'oxford'),
    p = oxdata(varargin{:});
end 
if strcmp(dset,'luna'),
    p = luna_data(varargin{:});
end 
if strcmp(dset,'yale'),
    [p,psall] = yale_data(varargin{:});
end 
if strcmp(dset,'flanker'),
    [p,psall] = flanker_data(varargin{:});
end 
if strcmp(dset,'sstxui'),
    [p,psall] = sstxui_data(varargin{:});
end 
if strcmp(dset,'stopgo'),
    p = stopgo_data(varargin{:});
end 
if strcmp(dset,'oddball'),
    p = oddball_data(varargin{:});
end 

if strcmp(dset,'artificial'),
    p = MDSadata(varargin{:});
end 
if ~isfield(p,'s'), 
    disp('WARNING! Dataset not found:');
    disp(dset);
    assert(false);
end 
 
names = {'rats','hw','open','bms','attention','spDCM','luna','yale'};
if ~isfield(p, 'C'), 
    M = size(p.s(1).y,1); J = size(p.s(1).v,1) ;
    p.C = repmat( -eye(M)/10, [1, 1, J]);    
end
[M,~,J] = size(p.C);

if ~isfield(p,'opts'), p.opts = struct(); end
if ~isfield(p.opts, 'PC'),
    for j=1:J, 
        I = eye(M); 
        dmin = -0.05; 
        if j > 1, dmin = 0; end
        p.opts.PC(j).opts.IV =  vec(I) * [-.5, dmin] + vec(1-I) * [-1, 1] * .8;
    end    
end

p = psdd(p);
for j=1:numel(psall),
    psall{j} = psdd(psall{j});    
end

end

function p = psdd(p),
%% standardize the dataset to simplify further analysis. 
for j=1:length(p.s)
    y = p.s(j).y;
    y = spm_detrend(y',3);
    y = zscore(y);
    p.s(j).y = y'; 
end
%% 
end
function addpathif(dir)
if isdir(dir)
    addpath(dir)
end
end