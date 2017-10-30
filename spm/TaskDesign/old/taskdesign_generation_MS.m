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
  fprintf('Processing subject: %s\n', subject);
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
      ListName = ListFile(1).name;
      OGName = strcat('tabfile_',ListName);
      ListName = ListName(1:end-4);
      ListName = strcat('tabfile_',ListName,'.tsv');
      disp(OutputName)
      unix(sprintf('/home/groups/menon/lab_shared/software/.gem/ruby/2.4.0/bin/eprime2tabfile -o %s %s', OutputName, ListFile(1).name));    
     % unix(sprintf('mv %s %s',OGName, ListName)) 
    end
    ConvertedListFile = dir('tabfile_*.txt');
    disp(ConvertedListFile)
  end
end

% call python script which does the actual parsing and creates the
% task_design.m files that are output into the
% data/imaging/participants/PID/visitX/sessionX/task/ dir
unix(sprintf('python /oak/stanford/groups/menon/scsnlscripts/spm/TaskDesign/taskdesign_generation.py '))

