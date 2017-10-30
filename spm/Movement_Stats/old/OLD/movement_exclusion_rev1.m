function movement_exclusion(Config_File)

CurrentDir = pwd;

disp('------------------------------------------------------------------');
fprintf('If you get error messages, send the error message and config file to tianwenc@stanford.edu\n');
disp('------------------------------------------------------------------');
Config_File = strtrim(Config_File);
if ~exist(Config_File,'file')
  fprintf('Error: Cannot find the configuration file ... \n');
  return;
end
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

%-Configurations
ServerPath     = paralist.ServerPath;
PrepFolder     = paralist.PreprocessedFolder;
SubjList       = paralist.SubjectList;
CondList       = paralist.SessionList;
VoxSize        = paralist.VoxelSize;
ScanToScanCrit = paralist.ScanToScanCrit;

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

%--------------------------------------------------------------------------
Subjects = ReadList(SubjList);
NumSubjs = length(Subjects);
Conditions = ReadList(CondList);
NumConds = length(Conditions);
NumRuns = NumSubjs*NumConds;

RunIndex = zeros(NumRuns, 4);
[C1, C2] = meshgrid(1:NumSubjs, 1:NumConds);
RunIndex(:,1) = C1(:);
RunIndex(:,2) = C2(:);
RunIndex(:,3) = 1:NumRuns;
RunIndex(:,4) = 1;

MvmntDir = cell(NumRuns, 1);
%-overall max range | sum of max range | overall max scan to scan movement |
%-max of sum of scan to scan movement | # scans > 0.5 voxel w.r.t. max overall scan
%-to scan movement
MvmntStats = zeros(NumRuns, 8);

RunCnt = 1;
for iSubj = 1:NumSubjs
  for iCond = 1:NumConds
    YearId = ['20', Subjects{iSubj}(1:2)];
    MvmntDir{RunCnt} = fullfile(ServerPath, YearId, Subjects{iSubj}, 'fmri', ...
      Conditions{iCond}, PrepFolder);
    MvmntFile = fullfile(MvmntDir{RunCnt}, 'rp_I.txt');
    ZippedMvmntFile = fullfile(MvmntDir{RunCnt}, 'rp_I.txt.gz');
    GSFile = fullfile(MvmntDir{RunCnt}, 'VolumRepair_GlobalSignal.txt');
    ZippedGSFile = fullfile(MvmntDir{RunCnt}, 'VolumRepair_GlobalSignal.txt.gz');
    
    if exist(ZippedMvmntFile, 'file') || exist(ZippedGSFile, 'file')
      unix(sprintf('gunzip -fq %s', fullfile(MvmntDir{RunCnt}, '*.txt.gz')));
    else
      if ~exist(MvmntFile, 'file') || ~exist(GSFile, 'file')
        fprintf('Cannot find movement file or global signal file: %s\n', Subjects{iSubj});
        RunIndex(RunCnt, 4) = 0;
      else
        %-Load rp_I.txt
        load(MvmntFile);
        TransMvmnt = rp_I(:, 1:3);
        RotMvmnt = 50.*rp_I(:, 4:6);
        TotalDisp = sqrt(sum(TransMvmnt.^2, 2));
        ScanToScanTrans = abs(diff(TransMvmnt));
        ScanToScanTotalDisp = sqrt(sum(ScanToScanTrans.^2, 2));
        MvmntRange = range(TransMvmnt);
        
        MvmntStats(RunCnt, 1) = MvmntRange(1);
        MvmntStats(RunCnt, 2) = MvmntRange(2);
        MvmntStats(RunCnt, 3) = MvmntRange(3);
        %MvmntStats(RunCnt, 4) = max(range(TransMvmnt));
        MvmntStats(RunCnt, 4) = max(TotalDisp);
        %MvmntStats(RunCnt, 5) = max(ScanToScanTrans(:));
        MvmntStats(RunCnt, 5) = max(ScanToScanTotalDisp);
        %MvmntStats(RunCnt, 5) = sum(max(ScanToScanTrans, [], 2) > (ScanToScanCrit*VoxSize));
        MvmntStats(RunCnt, 6) = sum(ScanToScanTotalDisp > (ScanToScanCrit*VoxSize));
        mvnout_idx = (find(ScanToScanTotalDisp > (ScanToScanCrit*VoxSize)))'+1;
        
        g = load(GSFile);
        gsigma = std(g);
        gmean = mean(g);
        mincount = 5*gmean/100;
        %z_thresh = max( z_thresh, mincount/gsigma );
        z_thresh = mincount/gsigma;        % Default value is PercentThresh.
        z_thresh = 0.1*round(z_thresh*10); % Round to nearest 0.1 Z-score value
        zscoreA = (g - mean(g))./std(g);  % in case Matlab zscore is not available
        glout_idx = (find(abs(zscoreA) > z_thresh))';
        
        MvmntStats(RunCnt, 7) = length(glout_idx);
        
        union_idx = unique([1; mvnout_idx(:); glout_idx(:)]);
        MvmntStats(RunCnt, 8) = length(union_idx)/length(g)*100;
        
      end
    end
    RunCnt = RunCnt + 1;
  end
end

FullRunIndex = find(RunIndex(:,4) ~= 0);

if ~isempty(FullRunIndex)
  fid = fopen('MovementSummaryStats.txt', 'w+');
  fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'TASK', 'Scan_ID', ...
    'Range x', 'Range y', 'Range z', 'Max Displacement', ...
    'Max Scan-to-Scan Displacement', 'Num Scans > 0.5 Voxel Displacement', ...
    'Num Scans 5% Global Signal', '% of Volumes Repaired');
  for i = 1:length(FullRunIndex)
    fprintf(fid, '%s\t%s\t', Conditions{RunIndex(FullRunIndex(i), 2)}, ...
      Subjects{RunIndex(FullRunIndex(i), 1)});
    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', MvmntStats(RunIndex(FullRunIndex(i), 3), 1), ...
      MvmntStats(RunIndex(FullRunIndex(i), 3), 2), MvmntStats(RunIndex(FullRunIndex(i), 3), 3), ...
      MvmntStats(RunIndex(FullRunIndex(i), 3), 4), MvmntStats(RunIndex(FullRunIndex(i), 3), 5), ...
      MvmntStats(RunIndex(FullRunIndex(i), 3), 6), MvmntStats(RunIndex(FullRunIndex(i), 3), 7), ...
      MvmntStats(RunIndex(FullRunIndex(i), 3), 8));
  end
  fclose(fid);
  
  if length(FullRunIndex) < NumRuns
    fid = fopen('MovementMissingInfo.txt', 'w+');
    fprintf(fid, '%s\t%s\t%s\n', 'TASK', 'Scan_ID', 'DataDir');
    MissSet = setdiff(1:NumRuns, FullRunIndex);
    for i = 1:length(MissSet)
      fprintf(fid, 's\t%s\t%s\n', Conditions{RunIndex(MissSet(i), 2)}, ...
        Subjects{RunIndex(MissSet(i), 1)}, MvmntDir{MissSet(i)});
    end
    fclose(fid);
    
  end
else
  disp('None of the runs has rp_I.txt or global signal file');
end

cd(CurrentDir);
disp('------------------------------------------------------------------');
fprintf('Analysis is done!\n');
fprintf('Please check: MovementMissingInfo.txt (if any) for subjects that do not have movement files\n');
fprintf('Please check: MovementSummaryStats.txt for summary stats\n');
disp('------------------------------------------------------------------');
end

