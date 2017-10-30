function[] = build_contrasts_FuncConn(session)

% THIS VERSION: builds T-contrasts for seeded functional connectivity
% analysis. Assumes you have *zero* task conditions; allows any number
% of additional regressors (ROI timeseries, global signal, etc), as
% long as only one is "of interest".
% Author: catie 2007-02-20

% ORIGINAL VERSION:
% Function which creates contrasts for the statistical analysis.
% Author: Adnan Majid 2006-08-17
% Contact: dsamajid@stanford.edu

%keyboard;
contrastNames = {};
contrastVecs = {};
contrastTypes = {};
  
%%task_design = session.multi{1};
load('task_design.mat');

%%numConds = length(names);
numRegs = length(reg_vec);
%  vecLength = 2*numConds + numRegs + 1; % time deriv
%%vecLength = numConds + numRegs + 1;
vecLength = numRegs + 1;

base = zeros(1, vecLength);

contrastNames{1} = 'positive covariate';
contrastNames{2} = 'negative covariate';
  
%base(2*numConds+1:2*numConds+numRegs) = reg_vec; % time deriv
%%base(numConds+1:numConds+numRegs) = reg_vec;
base(1:numRegs) = reg_vec;  

contrastVecs{1} = base;
contrastVecs{2} = -base;
  
numTContrasts = length(contrastNames);
  
os = fopen('contrasts.txt', 'w');
fprintf(os, 'Contrasts for (functional connectivity) Statistical Analysis\n');
fprintf(os, 'Session: %s\n', sess_name);
for i=1:length(contrastNames)

  type = 'T';

  fprintf(os, '\n%s%d:\t%s\n', type, i, contrastNames{i});
  for j=1:length(contrastVecs{i}(:,1))
    fprintf(os, '\t[');
    fprintf(os, ' %d', contrastVecs{i}(j,:));
    fprintf(os, ' ]\n');
  end
end

fclose(os);
save contrasts.mat contrastNames contrastVecs numTContrasts;
