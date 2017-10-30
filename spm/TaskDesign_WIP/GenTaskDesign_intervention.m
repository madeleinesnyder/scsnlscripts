function GenTaskDesign_intervention(config_file)
%-------------------------------------------------------------------------------------
% GenTaskDesign_intervention(GenTaskDesign_intervention_config.txt)
%
%These values are specified in the config file

% subj_list     : Name of the text file containing the list of subjects.
%                 Specify the full path where the text file is.  If only the file name
%                 is specified, it is assumed that file exists in one of
%                 the matlab search paths. This file needs to contain both
%                 scan ID and subject number.
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
%
%TASK_DURATION  : -1 if you want to use the reaction time as the task duration else 
%                   specify a constant value you want assigned as task
%                   duration.
%
%-------------------------------------------------------------------------------------
%
%       OUTPUT:  
%               1)This program puts a task design file (named as you
%               specified) in each subjects 'task design' folder.  This file
%               contains the onsets and durations of all 'correct experimental'
%               'correct control' 'incorrect experimental' and 'incorrect control'
%               trials 
%               2) This program also outputs a file in the folder you ran
%               the script from.  This file contains a list of what  task
%               design files the program created 
%




config_file=strtrim(config_file);

if ~exist(config_file, 'file')
    fprintf('Cannot find the configuration file\n');
    return;
end

% get the variables from the config file
config_file = config_file(1:end-2);
eval(config_file);
clear config_file;


processingFolder=pwd;
keepingTrackSubject=[];
keepingTrackScan=[];
keepingTrackTask=[];
moveON=0;

RT_THRESHOLD = 150;  % Any response below 150msec should not be recorded

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(exist(subj_list, 'file') == 0)
  disp(sprintf('File does not exist: %s\n', subj_list));
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the list of subjects and their ID's

[scanNumb, subNumb] = GetSubjectsList_intervention(subj_list);

if(isempty(scanNumb))
  disp('NO SUBJECTS TO PROCESS');
  return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of procedures : They just act as STATE MACHINE labels.
task_accurate   = 1;
control_accurate = 2;
task_inaccurate    = 3;
control_inaccurate  = 4;
    % in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sess_name = upper(regexprep(expr_name, '_', ' '));

names{1}  = 'task_accurate';
names{2}  = 'control_accurate';
names{3}  = 'task_inaccurate';
names{4}  = 'control_inaccurate';
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
[exprname_col, subj_col,session_col,Block_col,equationAcc_col, equationResp_col, equationCresp_col,equationRT_col, onset_col,expName_col] = GetColIndices_intervention(str, expr_name);

RT_col  = equationRT_col;
         % data starts from the third line
         str = fgetl(fid);  
while(str ~= -1)
  
  linecnt   = 0;

  onsets    = {[] [] [] []};
	durations = {[] [] [] []};
    
    
  %This assigns which rows are Experimental condition, control condition or
  %rest
    addSub=[1,4,9,11,12,13,19,20,22,23,27,29];
    equals=[2,3,5,7,10,14,16,17,18,25,26,28];
    rest=[6,8,15,21,24,30];
    
  inds     = [0, regexpi(str, '\t'), length(str) + 1]; 

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % THE DATA PROCESSING STARTS HERE
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         

         block = str(inds(Block_col)+1:inds(Block_col+1)-1);  
         counter=1;
      %This loop goes through the list of 30 stimuli and categorizes them based on response
      
while(str2double(block)<30)
  
    moveON=0;
    if counter~=1
                str = fgetl(fid); 
   
    else
        counter=counter+3;
       onset_ref = str2double(str(inds(onset_col)+1:inds(onset_col+1)-1))/1000;
    end
        
  
        inds     = [0, regexpi(str, '\t'), length(str) + 1];
        currsubj = [str(inds(exprname_col)+1:inds(exprname_col+1)-1) '-' ...
              str(inds(subj_col)+1:inds(subj_col+1)-1)]; 
        subjectNumb=[str(inds(subj_col)+1:inds(subj_col+1)-1)];

        currExp=[str(inds(expName_col)+1:inds(expName_col+1)-1)]; 
 
         block = str(inds(Block_col)+1:inds(Block_col+1)-1);
         RT   = str2double(str(inds(equationRT_col)+1:inds(equationRT_col+1)-1));
         eResponse  = str2double(str(inds(equationResp_col)+1:inds(equationResp_col+1)-1));
         eCResponse = str2double(str(inds(equationCresp_col)+1:inds(equationCresp_col+1)-1));
          curr_onset = str2double(str(inds(onset_col)+1:inds(onset_col+1)-1))/1000;
    
    % Any response below 150msec should not be recorded
    % and considered as inaccurate response
    if(RT > 0 && RT < RT_THRESHOLD)
      acc = 0;
    end
  
    if(RT ==0)
	RT=9500;
    end
    
        if ~isempty(find(addSub==str2double(block),1))
           if eResponse== eCResponse
               
               curr_proc=task_accurate;
           else
               curr_proc=task_inaccurate;
           end
           
       elseif  ~isempty(find(equals==str2double(block),1))
            if eResponse==eCResponse
               
               curr_proc=control_accurate;
           else
               curr_proc=control_inaccurate;
            end
        elseif ~isempty(find(rest==str2double(block),1))
            moveON=6;
            
       end
    
    
    
    
    % ---------------------------------------------------------------------
    % Update the onset and the duration times
    
    if moveON<2
        
    
    if TASK_DURATION==-1
       
    
        
            durations{curr_proc} = [durations{curr_proc}  (RT/1000)];
       
            onsets{curr_proc}    = [onsets{curr_proc}      curr_onset - onset_ref];
    
            linecnt    = linecnt + 1;
            prev_proc  = curr_proc;

    else
        
    
        
            durations{curr_proc} = [durations{curr_proc}  TASK_DURATION];
       
            onsets{curr_proc}    = [onsets{curr_proc}      curr_onset - onset_ref];
    
            linecnt    = linecnt + 1;
            prev_proc  = curr_proc;

    end
    end
    

end

   if TASK_DURATION==-1
       
       
            
    if(str ~= -1)   % End of reading the file
      inds = [0 regexpi(str, '\t'), length(str) + 1];
      % The next subject's ID
      nextsubj = [str(inds(exprname_col)+1:inds(exprname_col+1)-1) '-' ...
              str(inds(subj_col)+1:inds(subj_col+1)-1)];

      % If the current subject and the next subject are different,
      % update the tsv file, re-initialize and get ready for the next subject
      
      
      
    end  
       
       
   else
         
    if(str ~= -1)   % End of reading the file
      inds = [0 regexpi(str, '\t'), length(str) + 1];
      % The next subject's ID
      nextsubj = [str(inds(exprname_col)+1:inds(exprname_col+1)-1) '-' ...
              str(inds(subj_col)+1:inds(subj_col+1)-1)];

      % If the current subject and the next subject are different,
      % update the tsv file, re-initialize and get ready for the next subject
      
     
 
    end 
   end
   
    
 
    ind    = find(subNumb==str2double(subjectNumb));
    currentScanId=scanNumb{ind};
    
  keepingTrackSubject=[keepingTrackSubject; subjectNumb];
  keepingTrackScan=[keepingTrackScan; currentScanId];
  keepingTrackTask=[keepingTrackTask;currExp];
  
  
    
   
    
        % ------------------------------------------------------------------------------------------- %
  % ------------------ TASK DESING PATH & FILE ------------------ %

  % If there is not ".m" at the end add it.
  if(isempty(regexp(task_dsgn, '.m$', 'ONCE')))
    task_dsgn = [task_dsgn '.m'];
  end

  
  if(ind<1)
    taskfile = [currsubj '_' task_dsgn];
    disp(sprintf('Data found for a subject that is not in the list of subjects : %s', currsubj));
  else

    taskpath = sprintf('/mnt/musk1/20%s/%s/fmri/%s/task_design/', currentScanId(1:2), currentScanId, currExp);
    if(~exist(taskpath, 'dir'))
      mkdir(taskpath);
    end
    taskfile = sprintf('%s/%s',taskpath, task_dsgn);

  end
  % ------------------------------------------------------------------------------------------- %
  % ------------------ CREATE THE TASK DESIGN M FILE ------------------ %
  disp(sprintf('Updating task design file: %s',taskfile));
  fidtask = fopen(taskfile, 'w');
  
    fprintf(fidtask, sprintf('sess_name     = ''%s'';\n\n',   currExp));
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
         str = fgetl(fid);  

 
end % END OF THE OUTER WHILE LOOP


taskName='TaskDesignCreated_';


 taskfile = sprintf('%s%s',taskName, task_dsgn);
 taskfile= regexprep(taskfile, '.m','.txt');

     fidfinal=fopen(taskfile, 'a+');

 
fprintf(fidfinal,'%s\n', datestr(now));
fprintf(fidfinal, '%s\n', '---------------------------------------------------------------------------');


for x=1:size(keepingTrackSubject,1)
    fprintf(fidfinal,'%s\t%s\t%s\n', keepingTrackSubject(x,:), keepingTrackScan(x,:),keepingTrackTask(x,:));
end
fclose(fidfinal);

