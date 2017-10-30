function[] = build_contrasts(session)

% Function which creates contrasts for the statistical analysis.
% Author: Adnan Majid 2006-08-17
% Contact: dsamajid@stanford.edu

contrastNames = {};
contrastVecs = {};
contrastTypes = {};

rest_exists = 0; %just in case the user forgets to write it into
reg_vec = [];    %the task_design file.

icont = 1;

if length(session) == 1
  
  task_design = session.multi{1};
  load(task_design);
  
  numConds = length(names);
  numRegs = length(reg_vec);
  vecLength = 2*numConds + numRegs + 1;
  base = zeros(1, vecLength);

  if ismember(1, reg_vec) % If a covariate of interest exists
    contrastNames{1} = 'positive covariate';
    contrastNames{2} = 'negative covariate';
    
    base(2*numConds+1:2*numConds+numRegs) = reg_vec;
    contrastVecs{1} = base;
    contrastVecs{2} = -base;
  
  else  % Regular T-contrasts
    for i=1:numConds
      for j=i:numConds
	if (i ~= j)
	  contrastNames{icont} = [names{i}, '-', names{j}];
	  contrastVecs{icont} = base;
	  contrastVecs{icont}(2*i - 1) = 1;
	  contrastVecs{icont}(2*j - 1) = -1;
	  icont = icont + 1;
	
	  contrastNames{icont} = [names{j}, '-', names{i}];
	  contrastVecs{icont} = -contrastVecs{icont-1};
	  icont = icont + 1;
	end
      end
    end

    % T contrasts for rests
    if (rest_exists == 1)
      for i=1:numConds
	contrastNames{icont} = [names{i}, '-rest'];
	contrastVecs{icont} = base;
	contrastVecs{icont}(2*i - 1) = 1;
	icont = icont + 1;

	contrastNames{icont} = ['rest-', names{i}];
	contrastVecs{icont} = -contrastVecs{icont-1};
	icont = icont + 1;
      end
    end
  end
  
  numTContrasts = length(contrastNames);

  % F contrasts
  fcontrasts = diag(ones(vecLength,1));
  fcontrasts = fcontrasts(1:vecLength-1,:);

  contrastNames{icont} = 'effects of interest';
  contrastVecs{icont} = fcontrasts;
  icont = icont + 1;

  for i=1:numConds
    contrastNames{icont} = names{i};
    contrastVecs{icont} = fcontrasts(2*i-1:2*i,:);
    icont = icont + 1;
  end

  os = fopen('contrasts.txt', 'w');
  fprintf(os, 'Contrasts for Single Session Statistical Analysis\n');
  fprintf(os, 'Session: %s\n', sess_name);
  
elseif length(session) == 2

  task_design_sess1 = session(:,1).multi{1};
  task_design_sess2 = session(:,2).multi{1};
  
  load(task_design_sess1);
  names1 = names;
  rest1 = rest_exists;
  sess1_name = sess_name;

  load(task_design_sess2);
  names2 = names;
  rest2 = rest_exists;
  sess2_name = sess_name;

  numConds = length(names);
  numRegs = length(reg_vec);
  
% NEEDS TO BE FIXED; HARD CODE MOVEMENT REGRESSORS; reg_vec never set
if (exist(session(:,1).multi_reg{1}))
  numRegs = 6; 
  end
  
  vecLength = 2*(2*numConds + numRegs + 1);
  ofst = 2*numConds + numRegs; %vector offset between session cond
  base = zeros(1, vecLength);

  if ismember(1, reg_vec)
    contrastNames{1} = 'positive covariate';
    contrastNames{2} = 'negative covariate';
    
    base(2*numConds+1:2*numConds+numRegs) = reg_vec;
    base(2*numConds+1+ofst:2*numConds+numRegs+ofst) = reg_vec;
    contrastVecs{1} = base;
    contrastVecs{2} = -base;
  
  else  % Regular T-contrasts    
    for i=1:numConds
      for j=i:numConds
	if (i ~= j)
	  contrastNames{icont} = ['s1_', names1{i}, '-s1_', names1{j}];
	  contrastVecs{icont} = base;
	  contrastVecs{icont}(2*i - 1) = 1;
	  contrastVecs{icont}(2*j - 1) = -1;
	  icont = icont + 1;
	  
	  contrastNames{icont} = ['s1_', names1{j}, '-s1_', names1{i}];
	  contrastVecs{icont} = -contrastVecs{icont-1};
	  icont = icont + 1;
	  
	  contrastNames{icont} = ['s2_', names2{i}, '-s2_', names2{j}];
	  contrastVecs{icont} = base;
	  contrastVecs{icont}(2*i - 1 + ofst) = 1;
	  contrastVecs{icont}(2*j - 1 + ofst) = -1;
	  icont = icont + 1;
	  
	  contrastNames{icont} = ['s2_', names2{j}, '-s2_', names2{i}];
	  contrastVecs{icont} = -contrastVecs{icont-1};
	  icont = icont + 1;
	  
	  contPair1 = ['s1_', names1{i}, '+s2_', names2{i}];
	  contPair2 = ['s1_', names1{j}, '+s2_', names2{j}];
	  contrastNames{icont} = [contPair1, '_MINUS_', contPair2];
	  contrastVecs{icont} = contrastVecs{icont-4} + ...
	      contrastVecs{icont-2};
	  icont = icont + 1;
	  
	  contrastNames{icont} = [contPair2, '_MINUS_', contPair1];
	  contrastVecs{icont} = -contrastVecs{icont-1};
	  icont = icont + 1;
	end
      end
    end
    
    for i=1:numConds
      for j=i:numConds
	if (i ~= j)
	  contPair1 = ['s1_', names1{i}, '-s1_', names1{j}];
	  contPair2 = ['s2_', names2{i}, '-s2_', names2{j}];
	  contrastNames{icont} = [contPair1, '_MINUS_', contPair2];
	  contrastVecs{icont} = base;
	  contrastVecs{icont}(2*i - 1) = 1;
	  contrastVecs{icont}(2*j - 1) = -1;
	  contrastVecs{icont}(2*i - 1 + ofst) = -1;
	  contrastVecs{icont}(2*j - 1 + ofst) = 1;
	  icont = icont + 1;
	  
	  contrastNames{icont} = [contPair2, '_MINUS_', contPair1];
	  contrastVecs{icont} = -contrastVecs{icont-1};
	  icont = icont + 1;
	end
      end
    end
    
    % T contrasts for rests
    if (rest1 == 1)
      for i=1:numConds
	contrastNames{icont} = ['s1_', names1{i}, '-rest'];
	contrastVecs{icont} = base;
	contrastVecs{icont}(2*i - 1) = 1;
	icont = icont + 1;
	
	contrastNames{icont} = ['rest-s1_', names1{i}];
	contrastVecs{icont} = -contrastVecs{icont-1};
	icont = icont + 1;
      end
    end
    
    if (rest2 == 1) %changed from rest1==2
      for i=1:numConds
	contrastNames{icont} = ['s2_', names2{i}, '-rest'];
	contrastVecs{icont} = base;
	contrastVecs{icont}(2*i - 1 + ofst) = 1;
	icont = icont + 1;
	
	contrastNames{icont} = ['rest-s2_', names2{i}];
	contrastVecs{icont} = -contrastVecs{icont-1};
	icont = icont + 1;
      end
    end
    
    if (rest1 == 1 & rest2 == 1)
      for i=1:numConds
	contPair = ['s1_', names1{i}, '+s2_', names2{i}];
	
	contrastNames{icont} = [contPair, '_MINUS_rest'];
	contrastVecs{icont} = base;
	contrastVecs{icont}(2*i - 1) = 1;
	contrastVecs{icont}(2*i - 1 + ofst) = 1;
	icont = icont + 1;
	
	contrastNames{icont} = ['rest_MINUS_', contPair];
	contrastVecs{icont} = -contrastVecs{icont-1};
	icont = icont + 1;
      end
      
      for i=1:numConds
	contrastNames{icont} = ['s1_', names1{i}, '-s2_', names2{i}];
	contrastVecs{icont} = base;
	contrastVecs{icont}(2*i - 1) = 1;
	contrastVecs{icont}(2*i - 1 + ofst) = -1;
	icont = icont + 1;
	
	contrastNames{icont} = ['s2_', names2{i}, '-s1_', names1{i}];
	contrastVecs{icont} = -contrastVecs{icont-1};
	icont = icont + 1;
      end  
    end
  end
  
  numTContrasts = length(contrastNames);

  % F contrasts
  fcontrasts = diag(ones(vecLength,1)); 
  fcontrasts = fcontrasts(1:vecLength-2,:);

  contrastNames{icont} = 'effects of interest';
  contrastVecs{icont} = fcontrasts; 
  icont = icont + 1;

  for i=1:numConds
    contrastNames{icont} = ['s1_', names1{i}];
    contrastVecs{icont} = fcontrasts(2*i-1:2*i,:);
    icont = icont + 1;
  end

  for i=1:numConds
    contrastNames{icont} = ['s2_', names2{i}];
    contrastVecs{icont} = fcontrasts(2*i-1+ofst:2*i+ofst,:);
    icont = icont + 1;
  end
  
  os = fopen('contrasts.txt', 'w');
  fprintf(os, 'Contrasts for 2-Session Statistical Analysis\n');
  fprintf(os, 'Session 1: %s\n', sess1_name);
  fprintf(os, 'Session 2: %s\n', sess2_name);
  
end

for i=1:length(contrastNames)

  if i <= numTContrasts
    type = 'T';
  else
    type = 'F';
  end

  fprintf(os, '\n%s%d:\t%s\n', type, i, contrastNames{i});
  for j=1:length(contrastVecs{i}(:,1))
    fprintf(os, '\t[');
    fprintf(os, ' %d', contrastVecs{i}(j,:));
    fprintf(os, ' ]\n');
  end
end

fclose(os);

save contrasts.mat contrastNames contrastVecs numTContrasts;
