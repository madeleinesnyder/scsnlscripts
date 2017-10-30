% Preprocess fMRI data with a specified pipeline
%__________________________________________________________________________
%-SCSNL, Tianwen Chen, 2011-12-02

function preprocessfmri(ConfigFile)

currentdir = pwd;

disp('==================================================================');
fprintf('Current directory: %s\n', currentdir);
fprintf('Script: %s\n', which('preprocessfmri.m'));
fprintf('Configfile: %s\n', ConfigFile);
disp('***Send error messages to: tianwenc@stanford.edu***');
fprintf('\n');

if ~exist(fullfile(ConfigFile), 'file')
    fprintf('Error: cannot find the configuration file... \n')
    return;
end

config = load(ConfigFile);
clear ConfigFile;

subject_i          = config.subject;
subjectlist        = strtrim(config.subjectlist);
runlist            = strtrim(config.runlist);
inputimgprefix     = strtrim(config.inputimgprefix);
wholepipeline      = strtrim(config.pipeline);
pipeline           = wholepipeline(1:end-length(inputimgprefix));
SPGRsubjectlist    = strtrim(config.SPGRsubjectlist);
TR                 = config.tr_val;
data_dir           = strtrim(config.data_dir);
project_dir        = strtrim(config.project_dir);

disp('-------------- Contents of the Parameter List --------------------');
disp(config);
disp('==================================================================');

%==========================================================================
%-Hard-coded configurations
%-Must notify tianwenc@stanford.edu if you make any changes below
data_type          = 'nii';
output_folder      = 'smoothed_spm8';
smooth_width       = [6 6 6];
boundingboxdim     = [-90 -126 -72; 90 90 108];
template_path      = '/oak/stanford/groups/menon/scsnlscripts/spm/BatchTemplates';
SPGR_folder        = 'anatomical';
try
    SPGRfilename   = strtrim(config.SPGRfilename);
catch e
    SPGRfilename   = 'watershed_spgr';
end

%==========================================================================

clear config;

if any(~ismember(data_type, {'nii', 'img'}))
  disp('Error: wrong data type specified');
  return;
end

if any(smooth_width < 0)
  disp('Error: smoothing kernel width cannot be negative');
  return;
end

if ~exist(template_path, 'dir')
  disp('Error: template folder does not exist!');
  return;
end

if ~isfloat(TR)
  disp('Error: TR must be a numerical float');
  return;
end

if ismember('f', wholepipeline)
  flipflag = 1;
else
  flipflag = 0;
end

subjectlist       = csvread(subjectlist,1);
subject           = subjectlist(subject_i);
subject           = char(pad(string(subject),4,'left','0'));
visit             = num2str(subjectlist(subject_i,2));
session           = num2str(subjectlist(subject_i,3));

numsubj           = length(subject);
runs              = ReadList(runlist);
numrun            = length(runs);

% are we using this still? JN 7/5/17
if ~isempty(SPGRsubjectlist)
  SPGRsubjects    = ReadList(SPGRsubjectlist);
  numSPGRsubj     = length(SPGRsubjects);
else
  SPGRsubject     = subject;
  numSPGRsubj     = numsubj;
end

if numSPGRsubj ~= numsubj
  disp('Number of functional subjects is not equal to the number of SPGR subjects');
  return;
end

numtotalrun = numsubj*numrun;
errmsg = cell(numtotalrun, 1);
errmsgflag = zeros(numtotalrun, 1);
totalrun_dir = cell(numtotalrun, 1);

volrepairflag = zeros(numtotalrun, 1);
volrepairdir = cell(numtotalrun, 1);

pipelinefamily = {'swar', 'swavr', 'swgcar', 'swgcavr', ...
  'swfar', 'swfavr', 'swgcfar', 'swgcfavr', ...
  'swaor', 'swgcaor', 'swfaor', 'swgcfaor'};

% if any(~ismember(wholepipeline, pipelinefamily))
%   disp('Error: unrecognized entire pipeline to be implemented');
%   return;
% end

spm('defaults', 'fmri');
spm_jobman('initcfg');
delete(get(0, 'Children'));

runcnt = 0;
for isubj = 1:numsubj
  fprintf('Processing subject: %s\n', subject);

  SPGRdir = fullfile(data_dir,SPGRsubject,['visit',visit],['session',session],SPGR_folder);
  
  SPGRfile_file = '';
  
  if ismember('c', wholepipeline)
    unix(sprintf('gunzip -fq %s', fullfile(SPGRdir, [SPGRfilename, '*.gz'])));
    listfile_file = dir(fullfile(SPGRdir, [SPGRfilename, '*']));
    SPGRfile_file = fullfile(SPGRdir, listfile_file(1).name);
  end

  for irun = 1:numrun
    runcnt = runcnt + 1;
    errcnt = 1;
    fprintf('---> run: %s\n', runs{irun});

    totalrun_dir{runcnt} = fullfile(data_dir, subject,['visit',visit],['session',session], 'fmri', ...
      runs{irun});


    % Put tmp directories in scratch in case  temp files get stuck.

    tmp_dir = fullfile('/scratch/users',getenv('LOGNAME'), 'tmp_files');
    if ~exist(tmp_dir, 'dir')
      mkdir(tmp_dir);
    end
      
    temp_dir = fullfile(tmp_dir, [subject,['visit',visit],['session',session], ...
	runs{irun},'_', tempname,'_', wholepipeline]);
    
    unnorm_dir = fullfile(totalrun_dir{runcnt}, 'unnormalized');
    if isempty(inputimgprefix)
      if ~exist(temp_dir, 'dir')
        mkdir(temp_dir);
      else
        unix(sprintf('/bin/rm -rf %s', fullfile(temp_dir, '*')));
      end
      unix(sprintf('cp -af %s %s', fullfile(unnorm_dir, ['I*', data_type, '*']), ...
        temp_dir));
      unix(sprintf('gunzip -fq %s', fullfile(temp_dir, 'I*nii.gz')));
    end
    
    imaging_path = '/data/imaging/participants/';
   output_dir = fullfile(project_dir,imaging_path,subject,['visit',visit],['session',session],'fmri',...
       runs{irun}, output_folder);
     
    volrepairdir{runcnt} = temp_dir;

    if ~exist(output_dir, 'dir')
      mkdir(output_dir);
    end
    
    output_log = fullfile(output_dir, 'log');
    if ~exist(output_log, 'dir')
      mkdir(output_log);
    end

    if ~isempty(inputimgprefix)
      if ~exist(temp_dir, 'dir')
        errmsg{runcnt}{errcnt} = sprintf('Directory does not exist: %s\n', temp_dir);
        disp(errmsg{runcnt}{errcnt});
        errcnt = errcnt + 1;
        errmsgflag(runcnt) = 1;
        continue;
      end
      listfile_file = dir(fullfile(temp_dir, 'meanI*'));
      if isempty(listfile_file)
        errmsg{runcnt}{errcnt} = sprintf('Error: no meanI_anon* image found when inputimgprefix is not empty');
        disp(errmsg{runcnt}{errcnt});
        errcnt = errcnt + 1;
        errmsgflag(runcnt) = 1;
        continue;
      else
        meanimg_file = fullfile(temp_dir, listfile_file(1).name);
      end
    end

    prevprefix = inputimgprefix;
    nstep = length(pipeline);

    for cnt = 1:nstep

      p = pipeline(nstep-cnt+1);

      switch p
        case 'r'
          listfile_file = dir(fullfile(temp_dir, [prevprefix, 'I*.gz']));
          if ~isempty(listfile_file)
            unix(sprintf('gunzip -fq %s', fullfile(temp_dir, [prevprefix, 'I*.gz'])));
          else
            [inputimg_file, selecterr] = preprocessfmri_selectfiles(temp_dir, prevprefix, data_type);
            if selecterr == 1
              errmsg{runcnt}{errcnt} = sprintf('Error: no scans selected');
              disp(errmsg{runcnt}{errcnt});
              errcnt = errcnt + 1;
              errmsgflag(runcnt) = 1;
              break;
            end
            preprocessfmri_realign(wholepipeline, currentdir, template_path, inputimg_file, temp_dir)
            unix(sprintf('/bin/rm -rf %s', fullfile(temp_dir, '*.mat')));
          end

          listfile_file = dir(fullfile(output_dir, ['rp_', prevprefix, 'I*.txt*.gz']));
          if ~isempty(listfile_file)
            unix(sprintf('gunzip -fq %', fullfile(output_dir, ['rp_', prevprefix, 'I*.txt*.gz'])));
          else
            listfile_file = dir(fullfile(output_dir, ['rp_', prevprefix, 'I*.txt']));
            if isempty(listfile_file)
              unix(sprintf('cp -af %s %s', fullfile(temp_dir, ['rp_', prevprefix, 'I*.txt']), output_dir));
            end
          end

          listfile_file = dir(fullfile(temp_dir, ['mean', prevprefix, 'I*', data_type]));
          meanimg_file = fullfile(temp_dir, listfile_file(1).name);


          if strcmpi(data_type, 'img')
            p = spm_select('ExtFPList', temp_dir, ['^r', prevprefix, 'I.*\.img']);
          else
            p = fullfile(temp_dir, ['r', prevprefix, 'I_anon.nii']);
          end
          vy = spm_vol(p);
          numscan = length(vy);
          disp('calculating the global signals ...');
          fid = fopen(fullfile(output_dir, 'VolumRepair_GlobalSignal.txt'), 'w+');
          for iscan = 1:numscan
            fprintf(fid, '%.4f\n', spm_global(vy(iscan)));
          end
          fclose(fid);

        case 'v'
          volflag = preprocessfmri_VolRepair(temp_dir, data_type, prevprefix);
          volrepairflag(runcnt) = volflag;
          nifti3Dto4D(temp_dir, prevprefix);
          unix(sprintf('gunzip -fq %s', fullfile(temp_dir, ['v', prevprefix, 'I*.gz'])));

          if volflag == 1
            disp('Skipping Art_Global (v) step ...');
            break;
          else
            unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'art_deweighted.txt'), output_dir));
            %unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'ArtifactMask.nii'), output_log));
            unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'art_repaired.txt'), output_log));
            unix(sprintf('mv -f %s %s', fullfile(temp_dir, '*.jpg'), output_log));
          end

         case 'o'
          volflag = preprocessfmri_VolRepair_OVersion(temp_dir, data_type, prevprefix);
          volrepairflag(runcnt) = volflag;
          %nifti3Dto4D(temp_dir, prevprefix);
          unix(sprintf('mv -f %s %s', fullfile(temp_dir, ['v', prevprefix, 'I.nii.gz']), fullfile(temp_dir, ['o', prevprefix, 'I.nii.gz'])));
          unix(sprintf('gunzip -fq %s', fullfile(temp_dir, ['o', prevprefix, 'I*.gz'])));


          if volflag == 1
            disp('Skipping Art_Global (o) step ...');
            break;
          else
            unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'art_deweighted.txt'), fullfile(output_dir, 'art_deweighted_o.txt')));
            %unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'ArtifactMask.nii'), output_log));
            unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'art_repaired.txt'), fullfile(output_log, 'art_repaired_o.txt')));
            unix(sprintf('mv -f %s %s', fullfile(temp_dir, '*.jpg'), output_log));
          end


        case 'f'
          preprocessfmri_FlipZ(temp_dir, prevprefix);

        case 'a'
          [inputimg_file, selecterr] = preprocessfmri_selectfiles(temp_dir, prevprefix, data_type);
          if selecterr == 1
            errmsg{runcnt}{errcnt} = sprintf('Error: no scans selected');
            disp(errmsg{runcnt}{errcnt});
            errcnt = errcnt + 1;
            errmsgflag(runcnt) = 1;
            break;
          end
          preprocessfmri_slicetime(wholepipeline, template_path, inputimg_file, flipflag, temp_dir, TR);

        case 'c'
          [inputimg_file, selecterr] = preprocessfmri_selectfiles(temp_dir, prevprefix, data_type);
          if selecterr == 1
            errmsg{runcnt}{errcnt} = sprintf('Error: no scans selected');
            disp(errmsg{runcnt}{errcnt});
            errcnt = errcnt + 1;
            errmsgflag(runcnt) = 1;
            break;
          end
          preprocessfmri_coreg(wholepipeline, template_path, data_type, SPGRfile_file, meanimg_file, temp_dir, inputimg_file, prevprefix);
          break;

        case 'w'
          [inputimg_file, selecterr] = preprocessfmri_selectfiles(temp_dir, prevprefix, data_type);
          if selecterr == 1
            errmsg{runcnt}{errcnt} = sprintf('Error: no scans selected');
            disp(errmsg{runcnt}{errcnt});
            errcnt = errcnt + 1;
            errmsgflag(runcnt) = 1;
            break;
          end
          preprocessfmri_normalize(wholepipeline, currentdir, template_path, boundingboxdim, [pipeline, inputimgprefix], inputimg_file, meanimg_file, temp_dir, SPGRfile_file);

        case 'g'
          listfile_file = dir(fullfile(SPGRdir, 'seg', '*seg_sn.mat'));
          if isempty(listfile_file)
            errmsg{runcnt}{errcnt} = sprintf('Error: no segmentation has been done, use preprocessfmri_seg.m');
            disp(errmsg{runcnt}{errcnt});
            errcnt = errcnt + 1;
            errmsgflag(runcnt) = 1;
            break;
          else
            if strcmp(data_type, 'img')
              imglist_file = dir(fullfile(temp_dir, [prevprefix, 'I*.img']));
              hdrlist_file = dir(fullfile(temp_dir, [prevprefix, 'I*.hdr']));
              num_file = length(imglist_file);
              for i_file = 1:num_file
                unix(sprintf('cp -af %s %s', fullfile(temp_dir, imglist_file(i_file).name), ...
                  fullfile(temp_dir, ['g', imglist_file(i_file).name])));
                unix(sprintf('cp -af %s %s', fullfile(temp_dir, hdrlist_file(i_file).name), ...
                  fullfile(temp_dir, ['g', hdrlist_file(i_file).name])));
              end
            else
              listfile_file = dir(fullfile(temp_dir, [prevprefix, 'I.nii']));
              unix(sprintf('cp -af %s %s', fullfile(temp_dir, listfile_file(1).name), ...
                fullfile(temp_dir, ['g', listfile_file(1).name])));
            end
          end

        case 's'
          [inputimg_file, selecterr] = preprocessfmri_selectfiles(temp_dir, prevprefix, data_type);
          if selecterr == 1
            errmsg{runcnt}{errcnt} = sprintf('Error: no scans selected');
            disp(errmsg{runcnt}{errcnt});
            errcnt = errcnt + 1;
            errmsgflag(runcnt) = 1;
            break;
          end
          preprocessfmri_smooth(wholepipeline, template_path, inputimg_file, temp_dir, smooth_width);

      end
      prevprefix = [pipeline((nstep-cnt+1):nstep), inputimgprefix];
      disp('------------------------------------------------------------');
    end

    if strcmp(prevprefix(1), 's')
      for iinter = 2:length(prevprefix)
        interprefix = prevprefix(iinter:end);
        listfile_file = dir(fullfile(temp_dir, [interprefix, 'I*']));
        num_file = length(listfile_file);
        for iinter_file = 1:num_file
          unix(sprintf('/bin/rm -rf %s', fullfile(temp_dir, listfile_file(iinter_file).name)));
        end
      end
      unix(sprintf('/bin/rm -rf %s', fullfile(temp_dir, '*.mat*')));
      unix(sprintf('gzip -fq %s', fullfile(temp_dir, [prevprefix, 'I*'])));
      unix(sprintf('gzip -fq %s', fullfile(output_dir, 'mean*I*')));
      unix(sprintf('gzip -fq %s', fullfile(temp_dir, 'mean*I*')));
      unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'mean*I*'), output_dir));
      unix(sprintf('mv -f %s %s', fullfile(temp_dir, [prevprefix, 'I*']), output_dir));
      unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'log', '*.mat'), fullfile(output_dir, 'log')));
      unix(sprintf('mv -f %s %s', fullfile(temp_dir, 'log', '*.pdf'), fullfile(output_dir, 'log')));
      unix(sprintf('/bin/rm -rf %s', fullfile(output_dir,'*.gz.gz*')));
      listfile_file = dir(fullfile(output_dir, '*.mat*'));
      if ~isempty(listfile_file)
        unix(sprintf('/bin/rm -rf %s', fullfile(output_dir, '*.mat*')));
      end
      listfile_file = dir(fullfile(output_dir, '*.jpg*'));
      if ~isempty(listfile_file)
        unix(sprintf('/bin/rm -rf %s', fullfile(output_dir, '*.jpg*')));
      end
      unix(sprintf('/bin/rm -rf %s', temp_dir));
    end
  end
  if all(ismember('sc', [pipeline, inputimgprefix]))
    unix(sprintf('gzip -fq %s', SPGRfile_file));
  end
end

cd(currentdir);

disp('==================================================================');
if sum(errmsgflag) == 0
  if ~strcmp(prevprefix(1), 's') && ismember('c', wholepipeline)
    disp('Please check coregistration quality');
  else
    disp('Preprocessing finished');
  end
else
  c = fix(clock);
  err_file = sprintf('ErrMsg_preprocessfmri_%d_%d_%d_%d_%d_%d.txt', c);
  fprintf('Please check: %s\n', err_file)
  errindex = find(errmsgflag == 1);
  fid = fopen(err_file, 'w+');
  for i = 1:length(errindex)
    fprintf(fid, '%s\n', totalrun_dir{errindex(i)});
    for j = 1:length(errmsg{errindex(i)})
      fprintf(fid, '---> %s\n', errmsg{errindex(i)}{j});
    end
  end
  fclose(fid);
end

if ismember('v', pipeline)
  if sum(volrepairflag) > 0
    disp('Please check: VolumeRepair_Flagged_Subjects_runs.txt for flagged subject_runs');
    flagfid = fopen('VolumeRepair_Flagged_Subjects_runs.txt', 'w');
    volrepindex = find(volrepairflag == 1);
    for i = 1:length(volrepindex)
      fprintf(flagfid, '%s\n', volrepairdir{volrepindex(i)});
    end
    fclose(flagfid);
  end
end

delete(get(0, 'Children'));
clear all;
close all;
disp('==================================================================');

end
