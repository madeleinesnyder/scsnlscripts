
% This script computes Contrast vectors
%
% To run Contrast_Generator:
% start matlab: ml7spm8, in matlab prompt
% >> contrast_generator(config_file)

%Use this script and config file to define contrast vectors for each contrast you want
% IndividualStats to do. Use for single or multi-run analyses.
%
% Input:    Edit the necessary parameters set out in the
%           Contrast_Geneartor_config.m file
%
% Output:   1) contrast.mat file that IndividualStats can use for
%              multirun experiments. The file has three parts:
%              ContrastNames, ContrastVecs numTContrasts.
%           2) contrast.txt file containing the final contrast vectors
%
% Details: The contrast.mat file that the IndividualStats script uses
% has a complicated structure, particularly for multi-run experiments.
% This script lets the author define simple contrast vectors for only the
% desired between-condition contrasts. These contrast vectors are then
% translated into the large complicated vectors expected by spm.
% In this script, condition order is based on the order specified in the
% task design file. The script assumes between-run effects are not of
% interest, so it collapses over runs.
%
% If you do not use this script, IndividualStats will create a default set
% of contrasts (often too many and not ones you're interested in, e.g.
% incorrect trials).
%
% Suggestion for task design files: If you have many conditions of
% non-interests, like inaccurate trials, arrange them as later numbers in
% your task design file.
%

% _________________________________________________________________________
% 2009 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: Contrast_Generator.m  $
% -------------------------------------------------------------------------


function contrast_generator(config_file)
%% Begin Program

fprintf('Contrast Generator\n');
disp('==================================================================');

%% Importing Basic information and contrast definitions from Contrast_Generator_config.m
%Check if config file is there
%config_file='Contrast_Generator_config.m';

config_file=strtrim(config_file);

if ~exist(config_file, 'file')
    fprintf('Cannot find the configuration file\n');
    return;
end

% get the variables from the config file
config = load(config_file);
clear config_file;

% Load parameters

nCon 		= config.numcontrasts
nrun 		= config.numruns
comWithinrun    = config.comparewithin
movFill         = config.movementcorrection
contrastNames   = config.contrastnames
contrast        = nest_cell_array(config.contrast,'double')
disp(contrast)

nBeta           = 2;  % Number of betas per condition = 2 for SPM
runFillend 	= 1;  % Final fill = 1 per run for SPM
runLen 		= nCon*nBeta + movFill; % run vector length without end fill
disp(runLen)
numTContrasts 	= length(contrastNames); %Number of contrasts. Check this. If wrong, look at your list again.


%% Generate contrast vectors composed of zeros
fprintf('Generating contrast vectors composed of zeros...\n');
contrastVecs = {};

% Vector for each run is made of:
%   2 components for each condition (one for each beta)
%   6 for movement corrections
%   1 for end-of-run
% With multiple runs, like runs A and B, the order of the vector
% is: A.conditions A.movement B.conds B.movt A.end-of-run B.end-of-run
vecLen = nrun*(nCon*nBeta + movFill)+ (nrun*runFillend);
for i = 1:numTContrasts
    contrastVecs{i} = zeros(1,vecLen);
end

%% Identify positions along the contrast vectors that represent your conditions.
fprintf('Identifying positions allong the contrast vectors that represent your conditions...\n');

% Remember, each condition has two adjacent betas and you're only
% interested in the first for each condition.

% Define each run's first position in the contrast vectors
runStart = [1];
for j = 1:(nrun-1)
    runStart = [runStart runLen*j+1];
end

% conditionPositions = [runStart(1) runStart(1)+2 ... runStart(2)
% runStart(2)+2 ...etc]
conditionPositions = []; %set up a vector to hold the position information
runConPos = [];  % I ADDED THIS DOES IT SEEM CORRECT?  JT: Sure, seems fine.
for i = 1:nrun
    for j = 1:nCon %make the vector of positions for run i
        runConPos(j) = runStart(i) + 2*(j-1);
    end
    conditionPositions = [conditionPositions runConPos]; %concat with previous runs
end

% As a check length(conditionPositions) = nCon*nrun

fprintf('Placing the contrasts defined above into the contrast vectors of zero...\n');
disp('------------------------------------------------------------------');

if comWithinrun==1,
    %% Place the contrasts defined above into the contrast vectors of zeros
% Use this version if the run have the idential conditions
fprintf('Makivzvng contrasts vector for runs with identical conditions with no comparison between runs...\n');
    for k = 1:numTContrasts
        contrastVecs{k}(conditionPositions) = repmat(contrast{k},1,nrun);
    end
elseif comWithinrun==0;
    %% Place the contrasts defined above into the contrast vectors of zeros
% Use this version if the run are different and you want to compare
% between them
fprintf('Making contrasts vector for contrasts compairing conditions between runs...\n');
    for k = 1:numTContrasts
         contrastVecs{k} (conditionPositions) = contrast{k};
    end

% Check the results before outputing the values
else
    fprintf('comWithinrun variable has not been defined correctly... please check configuration file!');
    return;
end
%% Create .mat and .txt files
save(strcat(config.project_dir,'/data/imaging/participants/7605/visit3/session1/fmri/stats_spm8/contrasts.mat'),'contrastNames','contrastVecs','numTContrasts')

fid = fopen(strcat(config.project_dir,'/data/imaging/participants/7605/visit3/session1/fmri/stats_spm8/contrasts.txt'), 'wt');

for i = 1: numTContrasts
    fprintf(fid, '%d ', i);
    fprintf(fid, contrastNames{i} );
    fprintf(fid, '\n[ ');
    fprintf(fid,'%d ',contrastVecs{i});
    fprintf(fid, ']\n\n');

end
fclose(fid);
