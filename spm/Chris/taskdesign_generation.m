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

config=load(ConfigFile);

clear ConfigFile;
ServerPath         = strtrim(config.rawdir);
SubjectList        = strtrim(config.subject);
SessionList        = strtrim(config.sessionlist);
EDatFolder         = strtrim(config.edatfolder);
GlobalOnsetColName = strtrim(config.globalonsetcolumn);
GlobalOnsetShift   = config.globalonsetshift;
OnsetCols          = strtrim(config.onsetcolumn);
OffsetCols         = strtrim(config.offsetcolumn);
DurationCol        = strtrim(config.durationcolumn);
PresetDuration     = config.presetduration;
TypeColName        = strtrim(config.tasktypecolumn);
TaskNames          = strtrim(config.taskname);
MatchTaskName      = config.matchtaskname;
SessionName        = strtrim(config.sessionname);
TaskDesignFileName = strtrim(config.taskdesignfilename);
RestExist          = config.restexist;

disp('-------------- Contents of the Parameter List --------------------');
disp(config);
fprintf('\n');
clear config;

TypeColName = ReadList(TypeColName);
Subjects = ReadList(SubjectList);
NumSubj  = length(Subjects);
Sessions = ReadList(SessionList);
NumSess  = length(Sessions);
TotalSess = NumSubj*NumSess;
TaskNames = ReadList(TaskNames);
if isempty(SessionName)
  SessionName = Sessions;
else
  SessionName = ReadList(SessionName);
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

OKList = ones(TotalSess, 1);

DataDir = cell(TotalSess, 1);
TaskDesignDir = cell(TotalSess, 1);
ErrMsg = cell(TotalSess, 1);

SessCnt = 1;
for iSubj = 1:NumSubj
  disp('----------------------------------------------------------------');
  fprintf('Processing subject: %s\n', Subjects{iSubj});
  YearId = ['20', Subjects{iSubj}(1:2)];
  for iSess = 1:NumSess
    fprintf('---> Session: %s\n', Sessions{iSess});
    DataDir{SessCnt} = fullfile(ServerPath, YearId, Subjects{iSubj}, 'fmri', ...
      Sessions{iSess}, EDatFolder);
    TaskDesignDir{SessCnt} = fullfile(ServerPath, YearId, Subjects{iSubj}, 'fmri', ...
      Sessions{iSess}, 'task_design');
    cd(DataDir{SessCnt});
    ListFile = dir('*.txt');
    ConvertedListFile = dir('tabfile_*.txt');
    
    %----------------------------------------------------------------------
    %-Check the txt file for conversion
    if length(ListFile) - length(ConvertedListFile) > 1
      ErrMsg{SessCnt} = sprintf('Warning: more than 1 unconverted .txt file found in the Edat folder');
      disp(ErrMsg{SessCnt});
      OKList(SessCnt) = 0;
      SessCnt = SessCnt + 1;
      continue;
    end
    if isempty(ListFile)
      ErrMsg{SessCnt} = sprintf('Warning: no .txt file found in the Edat folder');
      disp(ErrMsg{SessCnt});
      OKList(SessCnt) = 0;
      SessCnt = SessCnt + 1;
      continue;
    end
    %----------------------------------------------------------------------
    
    %-Convert the txt file if not exist and read in the data
    if isempty(ConvertedListFile)
      OutputName = ['tabfile_', ListFile(1).name];
      unix(sprintf('eprime2tabfile -o %s %s', OutputName, ListFile(1).name));
    end
    ConvertedListFile = dir('tabfile_*.txt');
    
    fid = fopen(ConvertedListFile(1).name);
    RawData = {};
    cnt = 1;
    while ~feof(fid)
      linedata = textscan(fgetl(fid), '%s', 'Delimiter', '\t');
      RawData(cnt, :) = linedata{1}; %#ok<*AGROW>
      cnt = cnt+1;
    end
    fclose(fid);
    
    ColNames = RawData(1,:);
    
    %-Get global onset
    if ~isempty(GlobalOnsetColName)
      ColMatch = strcmpi(GlobalOnsetColName, ColNames);
      [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
      ErrMsg{SessCnt} = ErrMsgTmp;
      if ErrorFlag == 1
        disp(ErrMsg{SessCnt});
        OKList(SessCnt) = 0;
        SessCnt = SessCnt + 1;
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
        ErrMsg{SessCnt} = sprintf('Error: no onset or offset column found');
        disp(ErrMsg{SessCnt});
        OKList(SessCnt) = 0;
        SessCnt = SessCnt + 1;
        continue;
      end
      if length(OnsetColName) ~= length(OffsetColName)
        ErrMsg{SessCnt} = sprintf('Error: Number of onset column not equal to number of offset columns');
        disp(ErrMsg{SessCnt});
        OKList(SessCnt) = 0;
        SessCnt = SessCnt + 1;
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
        ErrMsg{SessCnt} = ErrMsgTmp;
        if ErrorFlag == 1
          disp(ErrMsg{SessCnt});
          OKList(SessCnt) = 0;
          SessCnt = SessCnt + 1;
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
      
      SessTaskNames = TaskNames;
      NumTask = length(SessTaskNames);
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
        ErrMsg{SessCnt} = sprintf(['Error: more than 1 onset or offset colums found ', ...
          'when the task type column is specified']);
        disp(ErrMsg{SessCnt});
        OKList(SessCnt) = 0;
        SessCnt = SessCnt + 1;
        continue;
      end
      
      OnsetColName = OnsetColName{1};
      ColMatch = strcmpi(OnsetColName, ColNames);
      [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
      ErrMsg{SessCnt} = ErrMsgTmp;
      if ErrorFlag == 1
        disp(ErrMsg{SessCnt});
        OKList(SessCnt) = 0;
        SessCnt = SessCnt + 1;
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
          ErrMsg{SessCnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{SessCnt});
            OKList(SessCnt) = 0;
            SessCnt = SessCnt + 1;
            continue;
          end
          OffsetColVec = str2double(RawData(2:end, ColMatch == 1));
          OffsetColVec = OffsetColVec(~isnan(OffsetColVec));
          OffsetColVec = OffsetColVec - GlobalOnset;
          Durations = OffsetColVec - OnsetColVec;
        else
          ColMatch = strcmpi(DurationCol{1}, ColNames);
          [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
          ErrMsg{SessCnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{SessCnt});
            OKList(SessCnt) = 0;
            SessCnt = SessCnt + 1;
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
        TaskDesign(iTask).Name = SessTaskNames{iTask};
      end
      
      %-The TaskTypeColumn is not specified
    else
      NumTask = length(OnsetColName);
      if isempty(TaskNames{1})
        ErrMsg{SessCnt} = sprintf('Error: task names are empty when tasktypecolumn is also empty');
        disp(ErrMsg{SessCnt});
        OKList(SessCnt) = 0;
        SessCnt = SessCnt + 1;
        continue;
      else
        if length(TaskNames) ~= NumTask
          ErrMsg{SessCnt} = sprintf('Error: number of task names not equal to number of onset columns when tasktypecolumn is empty');
          disp(ErrMsg{SessCnt});
          OKList(SessCnt) = 0;
          SessCnt = SessCnt + 1;
          continue;
        end
      end
      if isempty(DurationCol{1}) && isempty(PresetDuration)
        if length(OffsetColName) ~= length(OnsetColName)
          ErrMsg{SessCnt} = sprintf('Error: number of onset columns not equal to number of offset columns');
          disp(ErrMsg{SessCnt});
          OKList(SessCnt) = 0;
          SessCnt = SessCnt + 1;
          continue;
        end
      else
        NumDur = length(DurationCol);
        if NumDur == 1
          ColMatch = strcmpi(DurationCol{1}, ColNames);
          [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
          ErrMsg{SessCnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{SessCnt});
            OKList(SessCnt) = 0;
            SessCnt = SessCnt + 1;
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
            ErrMsg{SessCnt} = sprintf('Error: number of duration columns not equal to number of onset colums');
            disp(ErrMsg{SessCnt});
            OKList(SessCnt) = 0;
            SessCnt = SessCnt + 1;
            continue;
          else
            InnerLoopExit = 0;
            for iDur = 1:NumDur
              ColMatch = strcmpi(DurationCol{iDur}, ColNames);
              [ErrorFlag, ErrMsgTmp] = CheckColMatch(ColMatch);
              ErrMsg{SessCnt} = ErrMsgTmp;
              if ErrorFlag == 1
                disp(ErrMsg{SessCnt});
                OKList(SessCnt) = 0;
                SessCnt = SessCnt + 1;
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
        ErrMsg{SessCnt} = ErrMsgTmp;
        if ErrorFlag == 1
          disp(ErrMsg{SessCnt});
          OKList(SessCnt) = 0;
          SessCnt = SessCnt + 1;
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
          ErrMsg{SessCnt} = ErrMsgTmp;
          if ErrorFlag == 1
            disp(ErrMsg{SessCnt});
            OKList(SessCnt) = 0;
            SessCnt = SessCnt + 1;
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
    
    if ~exist(TaskDesignDir{SessCnt}, 'dir')
      mkdir(TaskDesignDir{SessCnt});
    end
    DesignFile = fullfile(TaskDesignDir{SessCnt}, TaskDesignFileName);
    
    %-Write out task_design_*.m
    fid = fopen(DesignFile, 'w+');
    fprintf(fid, 'sess_name = ''%s'';\n', SessionName{iSess});
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
    
    fprintf(fid, 'save task_design.mat sess_name names onsets durations rest_exists');
    
    fclose(fid);
    
    SessCnt = SessCnt + 1;
    clear RawData;
  end
end

disp('------------------------------------------------------------------');
fprintf('Done!\n');
cd(CurrentDir);
if sum(OKList) ~= TotalSess
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
