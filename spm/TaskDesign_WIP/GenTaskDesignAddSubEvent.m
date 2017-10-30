function GenTaskDesignAddSubEvent(args)
%-------------------------------------------------------------------------------------
% GenTaskDesignAddSubEvent(subj_list, expr_name, data_file, task_dsgn, rest_exists)
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
% expr_name   = 'addition_event_1';
% data_file   = 'event_merged.txt';
% task_dsgn   = 'task_design_generated.m';
% rest_exists = 1;

RT_THRESHOLD = 150;  % Any response below 150msec should not be recorded

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
PROC_COMPLEX_ACC   = 1;
PROC_COMPLEX_INACC = 2;
PROC_SIMPLE_ACC    = 3;
PROC_SIMPLE_INACC  = 4;

TASK_DURATION = 5.5;    % in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sess_name = upper(regexprep(expr_name, '_', ' '));

names{1}  = 'complex accurate';
names{2}  = 'complex inaccurate';
names{3}  = 'simple accurate';
names{4}  = 'simple inaccurate';
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
ctrlcresp_col  ctrlresp_col  ctrlRT_col  stimulus_col jitter_col onset_col] = GetColIndices(str, expr_name);

RT_col  = addRT_col;
acc_col = addacc_col;
str = fgetl(fid);           % data starts from the third line

while(str ~= -1)
  
  linecnt   = 0;

  onsets    = {[] [] [] []};
	durations = {[] [] [] []};
  
  jittersum = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % THE DATA PROCESSING STARTS HERE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  NEWSUBJECT = 0;
  
  inds     = [0, regexpi(str, '\t'), length(str) + 1];
  currsubj = [str(inds(exprname_col)+1:inds(exprname_col+1)-1) '-' ...
              str(inds(subj_col)+1:inds(subj_col+1)-1)];
  onset_ref = str2double(str(inds(onset_col)+1:inds(onset_col+1)-1))/1000;

  while(~NEWSUBJECT)
    inds = [0 regexpi(str, '\t'), length(str) + 1];
    stim = str(inds(stimulus_col)+1:inds(stimulus_col+1)-1);
    RT   = str2double(str(inds(RT_col)+1:inds(RT_col+1)-1));
    acc  = str2double(str(inds(acc_col)+1:inds(acc_col+1)-1));
    
    jitter    = str2double(str(inds(jitter_col)+1:inds(jitter_col+1)-1))/1000;
    
    curr_onset = str2double(str(inds(onset_col)+1:inds(onset_col+1)-1))/1000;
    
    % Any response below 150msec should not be recorded
    % and considered as inaccurate response
    if(RT > 0 && RT < RT_THRESHOLD)
      acc = 0;
    end
    
    % this would be 1 in the case of simple math.
    simple = ~isempty([regexp(stim, '1\s*[\+\-]\s*\d+\s*=\s*(\d+)'), regexp(stim, '\d+\s*[\+\-]\s*1\s*=\s*(\d+)')]);
    
    if(simple == 0 && acc == 0)
      curr_proc = PROC_COMPLEX_INACC;
    
    elseif(simple == 0 && acc == 1)
      curr_proc = PROC_COMPLEX_ACC;

    elseif(simple == 1 && acc == 0)
      curr_proc = PROC_SIMPLE_INACC;
    
    elseif(simple == 1 && acc == 1)
      curr_proc = PROC_SIMPLE_ACC;
    end
    % ---------------------------------------------------------------------
    % Update the onset and the duration times
    if(linecnt ~= 0)
      durations{prev_proc} = [durations{prev_proc}  (curr_onset - prev_onset) - prevjitter];
    end
    onsets{curr_proc}    = [onsets{curr_proc}      curr_onset - onset_ref];
    
    jittersum  = jittersum + jitter;
    linecnt    = linecnt + 1;
    prev_proc  = curr_proc;
    prev_onset = curr_onset;
    prevjitter = jitter;
    
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
        durations{prev_proc} = [durations{prev_proc}  TASK_DURATION];
        NEWSUBJECT = 1;
      end
    else
      durations{prev_proc} = [durations{prev_proc}  TASK_DURATION];
      NEWSUBJECT = 1;
    end
  end % end of while NEWSUBJECT
  
  
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
    taskpath = sprintf('/fs/plum2_share1/20%s/%s/fmri/%s/task_design/', subj(1:2), subj, expr_name);
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
    fprintf(fidtask, sprintf('onsets{%d}    = [%s];\n',   col, num2str(onsets{col}, '%.2f, ')));
    fprintf(fidtask, sprintf('durations{%d} = [%s];\n\n', col, num2str(durations{col}, '%.2f, ')));
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
