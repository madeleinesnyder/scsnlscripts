function effconn_multisess_ppi_volrep(Config_File)

% This script performs psychophysiological interactions (PPI) analysis
% Configuration file: multisess_ppi_config.m.template
% _________________________________________________________________________
% 2009-2014 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: effconn_multisess_ppi_volrep.m 2009-04-09 $
% Tianwen Chen, 04/09/10
% 04/27/10: automatic searching for effects of interest contrast
%           generate contrasts.mat
% 02/23/2012: volume repair option
% 01/2014: multi session
% -------------------------------------------------------------------------

scsnl_id = '$Id: effconn_multisess_ppi_volrep.m 2014-02-18 v3 $';

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('PPI analysis started at %d/%02d/%02d %02d:%02d:%02d \n',c);
fprintf('%s \n', scsnl_id);
disp('==================================================================');
fname = sprintf('effconn_ppi_volrep_TEST-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

%-Check existence of the configuration file
Config_File = strtrim(Config_File);

if ~exist(Config_File,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off;
  return;
end

Config_File = Config_File(1:end-2);

%-Read individual stats parameters
currentdir = pwd;
eval(Config_File);
clear Config_File;

%-Load parameters
server_path = strtrim(paralist.processed_dir);
subjectlist = strtrim(paralist.subjectlist_file);
prep_pipeline = strtrim(paralist.pipeline_type);
stats_folder = strtrim(paralist.statsfolder_file);
%exp_sesslist = strtrim(paralist.session_list);
%exp_sessions = ReadList(sessionlist_file); %if we want a .txt, which we dont...
exp_sessions = paralist.session_list; %expect a cell array here
numsess = length(exp_sessions);

contrast_weight = paralist.contrastweights_value;
contrast_name = strtrim(paralist.contrast_type);

addpath(genpath('/home/fmri/fmrihome/SPM/spm8_scripts'));
addpath(genpath('/home/fmri/fmrihome/SPM/spm8'));

if ismember('roicenter_value', fields(paralist)) && ismember('roiradius_value', fields(paralist))
    use_sphere_rois = 1;
    roi_center = paralist.roicenter_value;
    roi_name = ReadList(paralist.roiname_list);
    roi_radius = paralist.roiradius_value;
    num_roi = length(roi_name);
    if ~isa(roi_radius,'cell')
        roi_radius = {roi_radius};
    end
    if ~isa(roi_center,'cell')
        roi_center = {roi_center};
    end
    if ( num_roi ~= length(roi_radius) ) || ( num_roi ~= length(roi_center))
      error('Numbers of ROI centers, radii, and names must all agree');
    end
elseif ismember('roi_file', fields(paralist))
    use_sphere_rois = 0;
    roi_file = ReadList(paralist.roifile_list);
    roi_name = ReadList(paralist.roiname_list);
    if length(roi_name) ~= length(roi_file)
      error('number of ROI files not equal to number of ROI names');
    end
    num_roi = length(roi_file);
else
    error('Please specify either ''roi_center'' and ''roi_radius'' for spherical ROI(s) or ''roi_file_list'' for pre-defined ROI(s)');
end

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

%-Read in lists
parent_folder = '';
subjects = multi_readlist(subjectlist);
numsubj = length(subjects{1});

ppi_xY = struct();
roi_info_all = struct();

ppi_u = cell(1,num_roi);
subdir_name = cell(1,num_roi);
ppi_name = cell(1,num_roi);


if use_sphere_rois
    for iROI = 1:num_roi
        %-Read in ROI information
        %ppi_xY.Ic is determined by the script
        ppi_xY(iROI).xyz = roi_center{iROI}';
        ppi_xY(iROI).name = roi_name{iROI};

        ppi_xY(iROI).def = 'sphere';
        ppi_xY(iROI).spec = roi_radius{iROI};

        %-Read in contrast weight matrix
        ppi_u{iROI}(:,1) = contrast_weight(:,1);
        ppi_u{iROI}(:,2) = ones(size(contrast_weight,1),1);
        ppi_u{iROI}(:,3) = contrast_weight(:,2);
        ppi_name{iROI} = ['(', roi_name{iROI}, ')', '*', '(', contrast_name, ')'];

        %-Folder name under *_PPI
        subdir_name{iROI} = [roi_name{iROI}, '_', contrast_name];

        %-Save roi info in each ppi stats folder for future reference
        roi_info_all(iROI).name = ppi_xY.name;
        roi_info_all(iROI).center = ppi_xY.xyz;
        roi_info_all(iROI).shape = ppi_xY.def;
        roi_info_all(iROI).spec = ppi_xY.spec;
        roi_info_all(iROI).unit = 'mm';
    end
else %TODO: sort out how to build the roi's from a nifti file
    for iROI = 1:num_roi

        roi_vol = spm_vol(roi_file{iROI});
        [roi_vals,roi_xyz] = spm_read_vols(roi_vol);
        %A = find(roi_vals);
        Q = find(roi_vals(:)>0); % linear matrix indices with ones
        ppi_xY(iROI).xyz = mean(roi_xyz(:,Q), 2); % mean of the roi is a "center"
        %ppi_xY(iROI).def = 'roi'; %TODO: 'mask' ???

        ppi_xY(iROI).rej = {'cluster'};
        ppi_xY(iROI).name = roi_name{iROI};
        ppi_xY(iROI).def = 'mask';
        ppi_xY(iROI).spec = roi_vol;

        %-Read in contrast weight matrix
        ppi_u{iROI}(:,1) = contrast_weight(:,1);
        ppi_u{iROI}(:,2) = ones(size(contrast_weight,1),1);
        ppi_u{iROI}(:,3) = contrast_weight(:,2);
        ppi_name{iROI} = ['(', roi_name{iROI}, ')', '*', '(', contrast_name, ')'];

        %-Folder name under *_PPI
        subdir_name{iROI} = [roi_name{iROI}, '_', contrast_name];

        %-Save roi info in each ppi stats folder for future reference
        roi_info_all(iROI).name = ppi_xY(iROI).name;
        roi_info_all(iROI).center = ppi_xY(iROI).xyz; %TODO
        roi_info_all(iROI).shape = ppi_xY(iROI).def;
        roi_info_all(iROI).spec = ppi_xY(iROI).spec;
        roi_info_all(iROI).unit = 'mm'; %TODO
    end

end

addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/ArtRepair'));

for iROI = 1:num_roi
    roi_info = roi_info_all(iROI);
    %-Begin processing
    for subcnt = 1:numsubj %Subject level loop
      if isempty(parent_folder)
        secondpart = subjects{1}{subcnt}(1:2);
        if str2double(secondpart) > 96
          pfolder = ['19' secondpart];
        else
          pfolder = ['20' secondpart];
        end
      else
        pfolder = parent_folder{1}{1};
      end
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>');
        fprintf('Processing subject: %s \n', subjects{1}{subcnt});
        stats_dir = fullfile(server_path, pfolder, subjects{1}{subcnt}, ...
          'fmri', 'stats_spm8', stats_folder);

        ppi_dir = fullfile(server_path, pfolder, subjects{1}{subcnt}, ...
          'fmri', 'stats_spm8', [stats_folder, '_PPI'], subdir_name{iROI});

        ppi_tempstats_dir = fullfile(ppi_dir, 'temp_stats');

        if ~exist(ppi_dir, 'dir')
          fprintf('Creating directory: %s\n', ppi_dir);
          mkdir(ppi_dir);
        else
          rmdir(ppi_dir, 's');
          mkdir(ppi_dir);
        end

        if ~exist(ppi_tempstats_dir, 'dir')
          mkdir(ppi_tempstats_dir);
        end

        cd(stats_dir);
        spm_original = load('SPM.mat');

        % Unzip raw images
        numimages = length(spm_original.SPM.xY.VY);
        %spm_original.SPM.xY.VY -- [total # of frames x 1] struct array
        ImgName = spm_original.SPM.xY.VY(1).fname;
        ImgPath = fileparts(ImgName);

        scan_number = [1 spm_original.SPM.nscan(1:end-1)];
        start_number = cumsum(scan_number);

        imgname = spm_original.SPM.xY.VY(1).fname;
        [pathstr, fname, fext] = fileparts(imgname);

        %get the session numbers for the provided sessions
        exp_sessions_numbers = zeros(numsess,1);

        if iscell(exp_sessions) %we are given session names --> find them based on image filenames
            fname_sessions = {spm_original.SPM.xY.VY(start_number).fname}';
            for ix = 1:numsess
                pat = sprintf('/%s/',exp_sessions{ix});
                iExp = find(cell2mat(cellfun(@(x)length(regexpi(x,pat)),fname_sessions,'uni',0)));
                exp_sessions_numbers(ix) = iExp;
            end
        else %assume exp_sessions is an array with session indices
            exp_sessions_numbers = exp_sessions;
        end

        if strcmp(fext, '.nii')
          for i = 1:length(start_number) %unzip images from all runs in SPM.mat
            imgname = spm_original.SPM.xY.VY(start_number(i)).fname;
            unix(sprintf('gunzip -fq %s', [imgname, '.gz']));
          end
        else
          for i = 1:numimages
            imgname = spm_original.SPM.xY.VY(i).fname;
            unix(sprintf('gunzip -fq %s', [imgname, '.gz']));
          end
        end

        % Find the number of effects of interest contrast
        numcons = length(spm_original.SPM.xCon);
        count = 1;
        for i = 1:numcons
          %spm_original.SPM.xCon(i).name
          if strcmp(spm_original.SPM.xCon(i).name, 'effects of interest')
            ppi_xY(iROI).Ic = count;
            fprintf('contrast %d is the effects of interest\n', count);
            break;
          else
            count = count + 1;
          end
        end

        if count > numcons
          error('Cannot find the contrast: effects of interest');
        end

        %-Extract BOLD in the seed region
        disp('---> Extracting Seed BOLD ...');
        for iSess = 1:length(exp_sessions) %session-specific loop
          cd(stats_dir)
          session_number = exp_sessions_numbers(iSess);
          ppi_xY(iROI).Sess = session_number;

          xY = ppi_regions(ppi_dir, ppi_xY(iROI));

          disp('---> Generating PPI regressors');
          %-Generate PPI regressors
          PPI = effconn_peb_ppi(fullfile(pwd, 'SPM.mat'), 'ppi', xY, ...
            ppi_u{iROI}, ppi_name{iROI}, 1, ppi_dir);

          cd(ppi_tempstats_dir);
          if iSess == 1
            %-Copy batch_stats.mat to ppi folder
            unix(sprintf('cp -a %s .', fullfile(stats_dir, 'batch_stats.mat')));

            %-Update the batch file
            disp('---> Updating batch file');
            load batch_stats.mat;

            %-Save roi information
            save roi_info.mat roi_info;
          end

          tempsess = matlabbatch{1}.spm.stats.fmri_spec.sess(1);

          %-Number of images in a session and its starting number
          num_sess_images = spm_original.SPM.nscan(session_number);
          sess_start_image = start_number(session_number);
          sess_end_image = sess_start_image + num_sess_images -1;
          matlabbatch{1}.spm.stats.fmri_spec.dir{1} = pwd;
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess) = tempsess;
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).scans = {};
          if strcmp(fext, '.nii')
            for i = sess_start_image:sess_end_image
              ii = (i+1) - sess_start_image;
              imgfile = [spm_original.SPM.xY.VY(i).fname, ',', num2str(ii)];
              matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).scans{ii} = imgfile;
            end
          else
            for i = sess_start_image:sess_end_image
              ii = (i+1) - sess_start_image;
              imgfile = [spm_original.SPM.xY.VY(i).fname, ',1'];
              matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).scans{ii} = imgfile;
            end
          end

          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi = {''};
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).multi_reg = {''};
          matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

          %-PPI interaction term
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress(1).name = ...
            ['PPI_', ppi_name{iROI}];
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress(1).val = PPI.ppi;

          %-Seed timeseries
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress(2).name = ...
            ['SeedBOLD_', roi_name{iROI}];
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress(2).val = PPI.Y;

          %-Contrast of task designs
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress(3).name = ...
            ['Psych_', contrast_name];
          matlabbatch{1}.spm.stats.fmri_spec.sess(iSess).regress(3).val = PPI.P;

        end %end session-specific loop

        matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(pwd, 'SPM.mat');

        matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(pwd, 'SPM.mat');
        matlabbatch{3}.spm.stats.con.consess = {''};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'PPI-Interaction';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 0 0 0];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl'; %replicate over sessions

        contrastNames{1} = ['PPI_', ppi_name{iROI}];
        contrastVecs{1} = [1 0 0 0];
        numTContrasts = 1;
        save contrasts.mat contrastNames contrastVecs numTContrasts;

        spm('defaults', 'fmri');
        %-Delete old batch file and SPM.mat file
        unix('/bin/rm -rf batch_stats.mat');
        save batch_stats_ppi matlabbatch
        %-Initialize the batch system
        disp('---> Run PPI estimation');
        spm_jobman('initcfg');
        delete(get(0,'Children'));
        %-Run analysis
        spm_jobman('run', 'batch_stats_ppi.mat');

        %ImgName = spm_original.SPM.xY.VY(1).fname;
        %ImgPath = fileparts(ImgName);
        fname_sessions = {spm_original.SPM.xY.VY(start_number).fname}';
        ImgPath = cellfun(@fileparts,fname_sessions,'uni',0);
        ImgPath = ImgPath(exp_sessions_numbers);
        %ImgPath = cellstr(ImgPath);
        cd(ppi_dir);
        scsnl_art_redo(ppi_tempstats_dir, prep_pipeline, ppi_dir, ImgPath);
        cd(ppi_tempstats_dir);
        delete('SPM.mat');
        unix(sprintf('/bin/cp *.mat %s', ppi_dir));
        unix(sprintf('/bin/cp *.jpg %s', ppi_dir));
        cd(ppi_dir);
        unix(sprintf('/bin/rm -rf %s', ppi_tempstats_dir));

        if strcmp(fext, '.nii')
          for i = 1:length(start_number)
            imgname = spm_original.SPM.xY.VY(start_number(i)).fname;
            unix(sprintf('gzip -fq %s', imgname));
          end
        else
          for i = 1:numimages
            imgname = spm_original.SPM.xY.VY(i).fname;
            unix(sprintf('gzip -fq %s', imgname));
          end
        end
    end

end %end outer ROI loop

cd(currentdir);
disp('------------------------------------------------------------------');
fprintf('Changing back to the directory: %s \n', currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('PPI analysis finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');

diary off;
clear all;
close all;
end
