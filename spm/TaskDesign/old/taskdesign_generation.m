%-Generate task design from edat files
%-Use 'eprime2tabfile' command in optimus-ep
%__________________________________________________________________________
%-SCSNL, Tianwen Chen, 2011/10/17

function taskdesign_generation(ConfigFile)

disp('==================================================================');
disp('taskdesign_generation.m is running');
fprintf('Current directory is: %s\n', pwd);
fprintf('Config file is: %s\n', ConfigFile);
disp('------------------------------------------------------------------');
disp('Send error messages to tianwenc@stanford.edu');
disp('==================================================================');
fprintf('\n');

CurrentDir = pwd;

ConfigFile = strtrim(ConfigFile);

if ~exist(ConfigFile,'file')
  fprintf('Error: cnnot find the configuration file ... \n');
  return;
end

config = load(ConfigFile);
clear ConfigFile;

subject_i          = config.subject;
subjectlist        = strtrim(config.subjectlist);
data_dir           = strtrim(config.raw_dir);
runlist            = strtrim(config.runlist);
EDatFolder         = strtrim(config.edat_dir);
GlobalOnsetColName = strtrim(config.globalonsetcolumn);
GlobalOnsetShift   = config.globalonsetshift;
OnsetCols          = strtrim(config.onsetcolumn);
OffsetCols         = strtrim(config.offsetcolumn);
DurationCol        = strtrim(config.durationcolumn);
PresetDuration     = config.presetduration;
TypeColName        = strtrim(config.tasktypecolumn);
TaskNames          = strtrim(config.taskname);
MatchTaskName      = config.matchtaskname;
runname            = strtrim(config.runname);
TaskDesignFileName = strtrim(config.taskdesignfilename);
RestExist          = config.restexist;
project_dir        = strtrim(config.project_dir);

disp('-------------- Contents of the Parameter List --------------------');
disp(config);
fprintf('\n');
clear config;

TypeColName        = ReadList(TypeColName);
subjectlist        = csvread(subjectlist,1);
subject           = subjectlist(subject_i);
subject           = char(pad(string(subject),4,'left','0'));
visit             = num2str(subjectlist(subject_i,2));
session           = num2str(subjectlist(subject_i,3));

numsubj            = length(subject);
runs 		   = ReadList(runlist);
numrun 		   = length(runs);
numtotalrun 	   = numsubj*numrun;
TaskNames 	   = ReadList(TaskNames);

if isempty(runname)
  runname = runs;
else
  runname = ReadList(runname);
end

[~, filestr, extstr] = fileparts(TaskDesignFileName);
if isempty(extstr)
  extstr = '.m';
  TaskDesignFileName = [filestr, extstr];
else
  if ~strcmpi(extstr, '.m')
    disp('Error: wrong type for task design file (shoul be *.m file)');
    return;
  end
end

%Nest MatchTaskName cell array
MatchTaskName = nest_cell_array(MatchTaskName);

%-Convert MatchTaskName to a cellstr
if ~isempty(MatchTaskName)
  if ~iscell(MatchTaskName)
    NumMatch = length(MatchTaskName);
    MatchTask = cell(NumMatch, 1);
    for i = 1:NumMatch
      MatchTask{i} = num2str(MatchTaskName(i));
    end
  else
    NumCells = length(MatchTaskName);
    for i = 1:NumCells
      if iscell(MatchTaskName{i})
        NumInnerCells = length(MatchTaskName{i});
        for j = 1:NumInnerCells
          if isnumeric(MatchTaskName{i}{j})
            MatchTaskName{i}{j} = num2str(MatchTaskName{i}{j});
          end
        end
      end
    end
    MatchTask = MatchTaskName;
  end
else
  MatchTask = MatchTaskName;
end

OKList = ones(numtotalrun, 1);

DataDir = cell(numtotalrun, 1);
TaskDesignDir = cell(numtotalrun, 1);
ErrMsg = cell(numtotalrun, 1);

runcnt = 1;
for isubj = 1:numsubj
  disp('----------------------------------------------------------------');
  fprintf('Processing subject: %s\n', subject(isubj));
  for irun = 1:numrun
    fprintf('---> run: %s\n', runs{irun});
    DataDir{runcnt} = fullfile(data_dir, subject,['visit',visit],['session',session], 'fmri', ...
      runs{irun}, EDatFolder);
    TaskDesignDir{runcnt} = fullfile(project_dir,'data/imaging/participants', subject, ...
      ['visit',visit],['session',session],'fmri',runs{irun});
    cd(DataDir{runcnt});
    ListFile = dir('*.txt');
    ConvertedListFile = dir('tabfile_*.txt');
    
    %----------------------------------------------------------------------
    %-Check the txt file for conversion
    if length(ListFile) - length(ConvertedListFile) > 1
      ErrMsg{runcnt} = sprintf('Warning: more than 1 unconverted .txt file found in the Edat folder');
      disp(ErrMsg{runcnt});
      OKList(runcnt) = 0;
      runcnt = runcnt + 1;
      continue;
    end
    if isempty(ListFile)
      ErrMsg{runcnt} = sprintf('Warning: no .txt file found in the Edat folder');
      disp(ErrMsg{runcnt});
      OKList(runcnt) = 0;
      runcnt = runcnt + 1;
      continue;
    end
    %----------------------------------------------------------------------
    
    %-Convert the txt file if not exist and read in the data
    if isempty(ConvertedListFile)
      OutputName = ['tabfile_', ListFile(1).name];
      disp(OutputName)
      unix(sprintf('/home/groups/menon/lab_shared/software/.gem/ruby/2.4.0/bin/eprime2tabfile -o %s %s', OutputName, ListFile(1).name));
    end
    ConvertedListFile = dir('tabfile_*.txt');
    disp(ConvertedListFile)
    
    fid = fopen(ConvertedListFile(1).name);
    RawData = {};
    cnt = 1;
    while ~feof(fid)
      Linedata = tdfread(ConvertedListFile(1).name); % MS edits oct 18 2017 (testing, creates the 02ex garbage)
      linedata = textscan(fgetl(fid), '%s', 'Delimiter', '\t');
      RawData(cnt, :) = linedata{1}; %#ok<*AGROW> % This is a 1 row X column header of the textfile
      cnt = cnt+1;
    end
    fclose(fid);
    RawData = strrep(RawData,'""','');
    ColNames = RawData(1,:);
    
    %-Get global onset
    if ~isempty(GlobalOnsetColName)
      ColMatch = strcmpi(GlobalOnsetColName, ColNames);
      [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
      ErrMsg{runcnt} = ErrMsgTmp;
      if ErrorFlag == 1
        disp(ErrMsg{runcnt});
        OKList(runcnt) = 0;
        runcnt = runcnt + 1;
        continue;
      end
      GlobalOnset = str2double(RawData(2:end, ColMatch == 1));
    end
    
    GlobalOnset = GlobalOnset(~isnan(GlobalOnset));
    GlobalOnset = GlobalOnset(1) + GlobalOnsetShift;
    GlobalOnset= double(GlobalOnset);
    
    TaskDesign = struct();
    
    OnsetColName  = ReadList(OnsetCols);
    OffsetColName = ReadList(OffsetCols);
    DurationCol   = ReadList(DurationCol);
    
    if isempty(DurationCol{1}) && isempty(PresetDuration)
      if isempty(OnsetColName{1}) || isempty(OffsetColName{1})
        ErrMsg{runcnt} = sprintf('Error: no onset or offset column found');
        disp(ErrMsg{runcnt});
        OKList(runcnt) = 0;
        runcnt = runcnt + 1;
        continue;
      end
      if length(OnsetColName) ~= length(OffsetColName)
        ErrMsg{runcnt} = sprintf('Error: Number of onset column not equal to number of offset columns');
        disp(ErrMsg{runcnt});
        OKList(runcnt) = 0;
        runcnt = runcnt + 1;
        continue;
      end
    end
    
    %-Get task types if exists
    if ~isempty(TypeColName{1})
      ColMatch = zeros(size(RawData, 2), 1);
      NumTaskType = length(TypeColName);
%       if NumTaskType > 2
%         error('The maximum number of task columns is 2');
%       end
      TypeColVec = [];
      for iTaskType = 1:NumTaskType
        TempMatch = strcmpi(TypeColName{iTaskType}, ColNames);
        [ErrorFlag, ErrMsgTmp] = CheckColMatch(TempMatch);
        ErrMsg{runcnt} = ErrMsgTmp;
        if ErrorFlag == 1
          disp(ErrMsg{runcnt});
          OKList(runcnt) = 0;
          runcnt = runcnt + 1;
          continue;
        end
        ColMatch = ColMatch + TempMatch(:);
        TypeColVec = [TypeColVec, RawData(2:end, TempMatch)];
      end
      
      %-Deal with file type stimuli (discarding file paths and extensions)
      for iTaskType = 1:NumTaskType
        NumTotalType = length(TypeColVec(:, iTaskType));
        for iType = 1:NumTotalType
          [yn, loc] = ismember('\', TypeColVec{iType, iTaskType});
          if yn
            tmpname = TypeColVec{iType, iTaskType}(loc+1:end);
            [~, tmpname] = fileparts(tmpname);
            TypeColVec{iType, iTaskType} = tmpname;
          else
            [~, tmpname] = fileparts(TypeColVec{iType, iTaskType});
            TypeColVec{iType, iTaskType} = tmpname;
          end
        end
      end
      
      runTaskNames = TaskNames;
      NumTask = length(runTaskNames);
      TaskIndex = zeros(size(RawData, 1) - 1, NumTask);
      
      for iTask = 1:NumTask
        InnerMatchTask = MatchTask{iTask};
        if iscell(InnerMatchTask)
          if length(InnerMatchTask) ~= length(TypeColName)
            error('wrong specification of MatchTaskName');
          end
        else
          InnerMatchTask = cellstr(InnerMatchTask);
        end
        
        IntersectIndex = ones(size(RawData, 1) - 1, 1);
        for iCondTask = 1:size(TypeColVec, 2)
          if ~isempty(InnerMatchTask{iCondTask})
            TempIndx = strcmpi(InnerMatchTask{iCondTask}, TypeColVec(:, iCondTask));
            if isempty(TempIndx)
              error('empty trials');
            end
            IntersectIndex(~TempIndx) = 0;
          end
        end
        TaskIndex(logical(IntersectIndex), iTask) = 1;
        
      end
      
      if length(OnsetColName) > 1 || length(OffsetColName) > 1
        ErrMsg{runcnt} = sprintf(['Error: more than 1 onset or offset colums found ', ...
          'when the task type column is specified']);
        disp(ErrMsg{runcnt});
        OKList(runcnt) = 0;
        runcnt = runcnt + 1;
        continue;
      end
      
      OnsetColName = OnsetColName{1};
      ColMatch = strcmpi(OnsetColName, ColNames);
      [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
      ErrMsg{runcnt} = ErrMsgTmp;
      if ErrorFlag == 1
        disp(ErrMsg{runcnt});
        OKList(runcnt) = 0;
        runcnt = runcnt + 1;
        continue;
      end

      OnsetColVec = str2double(RawData(2:end, ColMatch == 1));
      OnsetColVec = OnsetColVec(~isnan(OnsetColVec));
      OnsetColVec = OnsetColVec - GlobalOnset;
      
      if ~isempty(PresetDuration)
        DurTmp = PresetDuration*1000;
      else
        if isempty(DurationCol{1})
          OffsetColName = OffsetColName{1};
          ColMatch = strcmpi(OffsetColName, ColNames);
          [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
          ErrMsg{runcnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{runcnt});
            OKList(runcnt) = 0;
            runcnt = runcnt + 1;
            continue;
          end
          OffsetColVec = str2double(RawData(2:end, ColMatch == 1));
          OffsetColVec = OffsetColVec(~isnan(OffsetColVec));
          OffsetColVec = OffsetColVec - GlobalOnset;
          Durations = OffsetColVec - OnsetColVec;
        else
          ColMatch = strcmpi(DurationCol{1}, ColNames);
          [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
          ErrMsg{runcnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{runcnt});
            OKList(runcnt) = 0;
            runcnt = runcnt + 1;
            continue;
          end
          Durations = str2double(RawData(2:end, ColMatch == 1));
          Durations = Durations(~isnan(Durations));
        end
      end
      
      TaskIndex = logical(TaskIndex);
      for iTask = 1:NumTask
        TaskDesign(iTask).Onset = OnsetColVec(TaskIndex(:, iTask))./1000;
        if isempty(PresetDuration)
          TaskDesign(iTask).Duration = Durations(TaskIndex(:, iTask))./1000;
        else
          TaskDesign(iTask).Duration = DurTmp./1000;
        end
        if isempty(TaskDesign(iTask).Duration)
          TaskDesign(iTask).Duration = 0;
        end
        if isempty(TaskDesign(iTask).Onset)
          TaskDesign(iTask).Onset = 1000;
          TaskDesign(iTask).Duration = 0;
        end
        TaskDesign(iTask).Name = runTaskNames{iTask};
      end
      
      %-The TaskTypeColumn is not specified
    else
      NumTask = length(OnsetColName);
      if isempty(TaskNames{1})
        ErrMsg{runcnt} = sprintf('Error: task names are empty when tasktypecolumn is also empty');
        disp(ErrMsg{runcnt});
        OKList(runcnt) = 0;
        runcnt = runcnt + 1;
        continue;
      else
        if length(TaskNames) ~= NumTask
          ErrMsg{runcnt} = sprintf('Error: number of task names not equal to number of onset columns when tasktypecolumn is empty');
          disp(ErrMsg{runcnt});
          OKList(runcnt) = 0;
          runcnt = runcnt + 1;
          continue;
        end
      end
      if isempty(DurationCol{1}) && isempty(PresetDuration)
        if length(OffsetColName) ~= length(OnsetColName)
          ErrMsg{runcnt} = sprintf('Error: number of onset columns not equal to number of offset columns');
          disp(ErrMsg{runcnt});
          OKList(runcnt) = 0;
          runcnt = runcnt + 1;
          continue;
        end
      else
        NumDur = length(DurationCol);
        if NumDur == 1
          ColMatch = strcmpi(DurationCol{1}, ColNames);
          [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
          ErrMsg{runcnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{runcnt});
            OKList(runcnt) = 0;
            runcnt = runcnt + 1;
            continue;
          end
          DurTmp = str2double(RawData(2:end, ColMatch == 1));
          DurTmp = DurTmp(~isnan(DurTmp));
          Durations = cell(NumTask, 1);
          for iDur = 1:NumTask
            Durations{iDur} = DurTmp;
          end
        else
          if NumDur ~= NumTask
            ErrMsg{runcnt} = sprintf('Error: number of duration columns not equal to number of onset colums');
            disp(ErrMsg{runcnt});
            OKList(runcnt) = 0;
            runcnt = runcnt + 1;
            continue;
          else
            InnerLoopExit = 0;
            for iDur = 1:NumDur
              ColMatch = strcmpi(DurationCol{iDur}, ColNames);
              [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
              ErrMsg{runcnt} = ErrMsgTmp;
              if ErrorFlag == 1
                disp(ErrMsg{runcnt});
                OKList(runcnt) = 0;
                runcnt = runcnt + 1;
                InnerLoopExit = 1;
                break;
              end
              DurTmp = str2double(RawData(2:end, ColMatch == 1));
              DurTmp = DurTmp(~isnan(DurTmp));
              Durations{iDur} = DurTmp;
            end
            if InnerLoopExit == 1
              continue;
            end
          end
        end
      end
      
      InnerLoopExit = 0;
      
      for iTask = 1:NumTask
        TaskOnsetColName = OnsetColName{iTask};
        ColMatch = strcmpi(TaskOnsetColName, ColNames);
        [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
        ErrMsg{runcnt} = ErrMsgTmp;
        if ErrorFlag == 1
          disp(ErrMsg{runcnt});
          OKList(runcnt) = 0;
          runcnt = runcnt + 1;
          InnerLoopExit = 1;
          break;
        end
        TaskOnsetColVec = str2double(RawData(2:end, ColMatch == 1));
        TaskOnsetColVec = TaskOnsetColVec(~isnan(TaskOnsetColVec));
        TaskOnsetColVec = TaskOnsetColVec - GlobalOnset;
        
        if isempty(DurationCol{1})
          TaskOffsetColName = OffsetColName{iTask};
          ColMatch = strcmpi(TaskOffsetColName, ColNames);
          [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
          ErrMsg{runcnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{runcnt});
            OKList(runcnt) = 0;
            runcnt = runcnt + 1;
            InnerLoopExit = 1;
            break;
          end
          TaskOffsetColVec = str2double(RawData(2:end, ColMatch == 1));
          TaskOffsetColVec = TaskOffsetColVec(~isnan(TaskOffsetColVec));
          TaskOffsetColVec = TaskOffsetColVec - GlobalOnset;
        end
        
        if ~isempty(PresetDuration)
          DurTmp = PresetDuration*1000;
        else
          if isempty(DurationCol{1})
            DurTmp = TaskOffsetColVec - TaskOnsetColVec;
          else
            DurTmp = Durations{iTask};
          end
        end
        
        TaskDesign(iTask).Onset = TaskOnsetColVec./1000;
        TaskDesign(iTask).Duration = DurTmp./1000;
        if isempty(TaskDesign(iTask).Duration)
          TaskDesign(iTask).Duration = 0;
        end
        if isempty(TaskDesign(iTask).Onset)
          TaskDesign(iTask).Onset = 1000;
          TaskDesign(iTask).Duration = 0;
        end
        TaskDesign(iTask).Name = TaskNames{iTask};
      end
      if InnerLoopExit == 1
        continue;
      end
    end
    
    if ~exist(TaskDesignDir{runcnt}, 'dir')
      mkdir(TaskDesignDir{runcnt});
    end
    DesignFile = fullfile(TaskDesignDir{runcnt}, TaskDesignFileName);
    disp(DesignFile)
    
    %-Write out task_design_*.m
    fid = fopen(DesignFile, 'w+');
    fprintf(fid, 'run_name = ''%s'';\n', runname{irun});
    fprintf(fid, '\n');
    for iTask = 1:NumTask
      fprintf(fid, 'names{%d} = ''%s'';\n', iTask, TaskDesign(iTask).Name);
      
      fprintf(fid, 'onsets{%d} = [ ', iTask);
      fprintf(fid, '%.2f ', TaskDesign(iTask).Onset(:)');
      fprintf(fid, '];\n');
      
      fprintf(fid, 'durations{%d} = [ ', iTask);
      fprintf(fid, '%.2f ', TaskDesign(iTask).Duration(:)');
      fprintf(fid, '];\n');
      
      fprintf(fid, '\n');
    end
    
    fprintf(fid, 'rest_exists = %d;\n', RestExist);
    fprintf(fid, '\n');
    
    fprintf(fid, 'save task_design.mat run_name names onsets durations rest_exists');
    
    fclose(fid);
    
    runcnt = runcnt + 1;
    clear RawData;
  end
end

disp('------------------------------------------------------------------');
fprintf('Done!\n');
cd(CurrentDir);
if sum(OKList) ~= numtotalrun
  BadIndex = find(OKList == 0);
  NumBad = length(BadIndex);
  fid = fopen('Generation_Errors.txt', 'w+');
  for i = 1:NumBad
    fprintf(fid, '%s\n', DataDir{BadIndex(i)});
    fprintf(fid, '%s\n', ErrMsg{BadIndex(i)});
    fprintf(fid, '\n');
  end
  fclose(fid);
  fprintf('***Please check the Generation_Errors.txt for unsuccessful generation!***\n');
else
  fprintf('All generations are successful!\n');
end
disp('------------------------------------------------------------------');
end

function [ErrorFlag, ErrMsg] = CheckColMatch(ColMatch)
ErrMsg = '';
ErrorFlag = 0;
if sum(ColMatch) == 0
  ErrMsg = sprintf('Error: no matched column found');
  ErrorFlag = 1;
end

if sum(ColMatch) > 1
  ErrMsg = sprintf('Error: more than 1 matched column found');
  ErrorFlag = 1;
end


end
