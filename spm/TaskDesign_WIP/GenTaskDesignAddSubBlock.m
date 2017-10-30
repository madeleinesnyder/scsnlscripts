function GenTaskDesignAddSubBlock(args)
%-------------------------------------------------------------------------------------
%Modified on April 8th 2010, to properly sort subtraction simple trials 
%
%GenTaskDesignAddSubBlock(subj_list, expr_name, data_file, task_dsgn, rest_exists)
%
% subj_list     : Name of the text file containing the list of subjects.
%                 Specify the full path where the text file is.  If only the file name
%                 is specified, it is assumed that file exists in one of
%                 the matlab search paths.
%  
% expr_name     : Session folder in p21.  For example if your data is in
%                '/fs/plum2_share1/2006/06-05-31.1/fmri/addition_block' 
%                 then 'addition_block' would be the the argument.
%
% data_file     : Text file exported from the e-merge file.
%
% task_dsgn     : Task design file name.
%
% rest_exists   : 0 if rest does not exists and 1 if rest exists.
%-------------------------------------------------------------------------------------

% clear;

% subj_list   = '/mnt/plum3_share2/Users/kumar/behavioral/subjectslist.txt';
% expr_name   = 'addition_block';
% data_file   = '/mnt/plum3_share2/Users/kumar/behavioral/test.txt';
% task_dsgn   = 'task_design_generated.m';
% rest_exists = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subj_list   = args{1};
expr_name   = args{2};
data_file   = args{3};
task_dsgn   = args{4};

rest_exists = 1;  % Default value of 1 ==> rest exists
nargs = size(args, 2);
if(nargs == 5)
  rest_exists = args{5};
end

if(exist(subj_list, 'file') == 0)
  disp(sprintf('File does not exist: %s\n', subj_list));
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the list of subjects and their ID's

[subjects, subjID] = GetSubjectID(subj_list, expr_name);

if(isempty(subjects))
  disp('NO SUBJECTS TO PROCESS');
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of procedures : They just act as STATE MACHINE labels.
PROC_COMPLEX = 1;
PROC_SIMPLE  = 2;
PROC_FIND    = 3;
PROC_FIXATE  = 4;

TASK_DURATION = 5.5;    % in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmpi(expr_name, 'addition_block'))
  sess_name = 'Addition Block';
elseif(strcmpi(expr_name, 'subtraction_block'))
  sess_name = 'Subtraction Block';
else
  disp('Invalid experiment name');
  return;
end
  
names{1}  = ['complex'];
names{2}  = ['simple'];
names{3}  = ['find'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(exist(data_file, 'file') == 0)
  disp(sprintf('Behavioral data file does not exist: %s\n', data_file));
  return;
end

fid = fopen(data_file);
str = '';
while(isempty(regexpi(str, 'Stimulus')))
  str = fgetl(fid);           % Read till the column names show up
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[exprname_col, subj_col, addacc_col addcresp_col addresp_col addRT_col  ctrlacc_col  ...
ctrlcresp_col  ctrlresp_col  ctrlRT_col  stimulus_col jitter_col onset_col] = GetColIndices(str, expr_name); %#ok<NASGU>

str = fgetl(fid);           % data starts from the third line

while(str ~= -1)
  
  prev_proc = 0;
  linecnt   = 0;
  blkcnt    = 0;

  onsets    = {[] [] []};
	durations = {[] [] []};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % THE DATA PROCESSING STARTS HERE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  NEWSUBJECT = 0;
  
  inds     = [0, regexpi(str, '\t'), length(str) + 1];
  currsubj = [str(inds(exprname_col)+1:inds(exprname_col+1)-1) '-' ...
              str(inds(subj_col)+1:inds(subj_col+1)-1)];

  if(strcmpi(currsubj, 'addblock_060414_jma-9'))
    disp('hai');
  end

  while(~NEWSUBJECT)
    inds = [0 regexpi(str, '\t'), length(str) + 1];
    stim = str(inds(stimulus_col)+1:inds(stimulus_col+1)-1);

    % classify the stimulus into one of the four tasks
    if(regexp(stim, '^\s*\*\s*$'));
        curr_proc = PROC_FIXATE;        
    else

      if(isempty(regexp(stim, '(\d+)\s*[\+\-]\s*(\d)\s*', 'ONCE')))
          curr_proc = PROC_FIND;
      else
          
        if(strcmpi(expr_name, 'addition_block'))
          if(~isempty([regexp(stim, '1\s*[\+\-]\s*\d+\s*=\s*(\d+)'), regexp(stim, '\d+\s*[\+\-]\s*1\s*=\s*(\d+)')]))
          curr_proc = PROC_SIMPLE;
        else
          curr_proc = PROC_COMPLEX;
        end
        else
           if(regexp(stim, '\d+\s*[\+\-]\s*1\s*=\s*(\d+)'));
          curr_proc = PROC_SIMPLE;
        else
          curr_proc = PROC_COMPLEX;
           end
        end
        
      end % END of checking for PROC_FIND
      
    end % END of checking for PROC_FIXATE
    
    
    % Update the onset and the duration times
    if (prev_proc ~= curr_proc)
        % skip this the first time
        if((linecnt ~= 0) && (prev_proc ~= PROC_FIXATE))
            durations{prev_proc} = [durations{prev_proc}  blkcnt  * TASK_DURATION];
        end

        if(curr_proc ~= PROC_FIXATE)
            onsets{curr_proc}    = [onsets{curr_proc}     linecnt * TASK_DURATION];
        end
        blkcnt = 0;
    end
    linecnt   = linecnt + 1;
    blkcnt    = blkcnt  + 1;
    prev_proc = curr_proc;
    
    % Read the next line
    str = fgetl(fid);
    
    if(str ~= -1)   % End of reading the file
      inds = [0 regexpi(str, '\t'), length(str) + 1];
      % The next subject's ID
      nextsubj = [str(inds(exprname_col)+1:inds(exprname_col+1)-1) '-' ...
              str(inds(subj_col)+1:inds(subj_col+1)-1)];

      % If the current subject and the next subject are different,
      % update the tsv file, re-initialize and get ready for the next subject
      if(~strcmpi(currsubj, nextsubj))
        durations{prev_proc} = [durations{prev_proc}  blkcnt  * TASK_DURATION];
        NEWSUBJECT = 1;
      end
    else
      durations{prev_proc} = [durations{prev_proc}  blkcnt  * TASK_DURATION];
      NEWSUBJECT = 1;
    end
  end % end of while
  
  
  %%%%%%%%%%%%%%%%%%% UPDATE THE TASK DESIGN FILE %%%%%%%%%%%%%%%%
  
  % ------------------------------------------------------------------------------------------- %
  % ------------------ TASK DESING PATH & FILE ------------------ %

  % If there is not ".m" at the end add it.
  if(isempty(regexp(task_dsgn, '.m$', 'ONCE')))
    task_dsgn = [task_dsgn '.m'];
  end

  subjind = find(ismember(subjID, currsubj));
  
  if(isempty(subjind))
    taskfile = [currsubj '_' task_dsgn];
    disp(sprintf('Data found for a subject that is not in the list of subjects : %s', currsubj));
  else
    subj     = char(subjects{subjind});
    taskpath = sprintf('/fs/musk1/20%s/%s/fmri/%s/task_design/', subj(1:2), subj, expr_name);
    if(~exist(taskpath, 'dir'))
      mkdir(taskpath);
    end
    taskfile = sprintf('%s/%s',taskpath, task_dsgn);
    subjects(subjind) = [];
    subjID  (subjind) = [];
  end

  % ------------------------------------------------------------------------------------------- %
  % ------------------ CREATE THE TASK DESIGN M FILE ------------------ %
  disp(sprintf('Updating task design file: %s',taskfile));
  fidtask = fopen(taskfile, 'w');
  
    fprintf(fidtask, sprintf('sess_name     = ''%s'';\n\n',   sess_name));
  for col = 1:size(onsets, 2)
    
    % If any of the conditions do not have a trial, set the onsets and
    % durations to be 0
    if(isempty(onsets{col}))
      onsets{col}    = [0];
      durations{col} = [0];
    end

    fprintf(fidtask, sprintf('names{%d}     = [''%s''];\n',   col, names{col}'));
    fprintf(fidtask, sprintf('onsets{%d}    = [%s];\n',   col, num2str(onsets{col}, '%.1f, ')));
    fprintf(fidtask, sprintf('durations{%d} = [%s];\n\n', col, num2str(durations{col}, '%.1f, ')));
  end

    fprintf(fidtask, sprintf('rest_exists  = %d;\n\n',   rest_exists));
  fprintf(fidtask, 'save task_design.mat sess_name names onsets durations rest_exists');
  fclose(fidtask);

  % ------------------------------------------------------------------------------------------- %

end % END OF THE OUTER WHILE LOOP
fclose(fid);  % END OF READING ONE SUBJECT'S INFORMATION


% Didn't found data for these subjects
for cnt = 1:length(subjects)
  disp(sprintf('Task design not computed for the subject : %s',subjects{cnt}));
end
