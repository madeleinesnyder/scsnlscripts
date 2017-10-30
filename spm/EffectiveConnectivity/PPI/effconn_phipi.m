function effconn_phipi (Config_File)

% This script performs physiophysiological interactions (PhiPI) analysis
% Configuration file: effconn_phipi_config.m.template
% _________________________________________________________________________
% 2009-2011 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: effconn_phipi.m 2011-07-08 $
% Tianwen Chen, 07/08/2011
% -------------------------------------------------------------------------

scsnl_id = '$ Id: effconn_phipi.m 2011-07-08 v1 $';

warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
disp('==================================================================');
fprintf('PhiPI analysis started at %d/%02d/%02d %02d:%02d:%02d \n',c);
fprintf('%s \n', scsnl_id);
disp('==================================================================');
fname = sprintf('effconn_phipi-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
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
parent_folder = strtrim(paralist.parent_dir);
subjectlist = strtrim(paralist.subjectlist_file);
stats_folder = strtrim(paralist.stats_dir);
roi_center = paralist.roicenter_list;
roi_name = strtrim(paralist.roiname_list);
roi_radius = paralist.roiradius_list;
session_number = paralist.session_value;
contrast_weight = paralist.contrastweights_list;
contrast_name = strtrim(paralist.contrast_type);

disp('-------------- Contents of the Parameter List --------------------');
disp(paralist);
disp('------------------------------------------------------------------');
clear paralist;

%-Read in lists
parent_folder = multi_readlist(parent_folder);
subjects = multi_readlist(subjectlist);
numsubj = length(subjects{1});

%-Read in ROI information
%ppi_xY.Ic is determined by the script
%-First ROI
ppi_xY.xyz = roi_center{1}(:);
ppi_xY.name = roi_name{1};
ppi_xY.Sess = session_number;
ppi_xY.def = 'sphere';
ppi_xY.spec = roi_radius(1);

%-Read in contrast weight matrix
ppi_u(:,1) = contrast_weight(:,1);
ppi_u(:,2) = ones(size(contrast_weight,1),1);
ppi_u(:,3) = contrast_weight(:,2);
ppi_name = ['(', roi_name{1}, ')', '*', '(', roi_name{2}, ')', '*', '(', contrast_name, ')'];

%-Folder name under *_PPI
subdir_name = [roi_name{1}, '_', roi_name{2}, '_', contrast_name];

%-Save roi info in each ppi stats folder for future reference
roi_info = struct();
roi_info(1).name = ppi_xY.name;
roi_info(1).center = ppi_xY.xyz;
roi_info(1).shape = ppi_xY.def;
roi_info(1).spec = ppi_xY.spec;
roi_info(1).unit = 'mm';

%-Begin processing
for subcnt = 1:numsubj
  if isempty(parent_folder{1}{1})
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
      'fmri', 'stats_spm8', [stats_folder, '_PPI'], subdir_name);
    if ~exist(ppi_dir, 'dir')
      fprintf('Creating directory: %s\n', ppi_dir);
      mkdir(ppi_dir);
    else
      rmdir(ppi_dir, 's');
      mkdir(ppi_dir);
    end
    cd(stats_dir);
    spm_original = load('SPM.mat');
    % Unzip raw images
    numimages = length(spm_original.SPM.xY.VY);
    scan_number = [1 spm_original.SPM.nscan(1:end-1)];
    start_number = cumsum(scan_number);
    imgname = spm_original.SPM.xY.VY(1).fname;
    [pathstr, fname, fext] = fileparts(imgname);
    if strcmp(fext, '.nii')
      for i = 1:length(start_number)
        imgname = spm_original.SPM.xY.VY(start_number(i)).fname;
        fexist = dir(imgname);
        if isempty(fexist)
          unix(sprintf('gunzip -fq %s', [imgname, '.gz']));
        end
      end
    else
      for i = 1:numimages
        imgname = spm_original.SPM.xY.VY(i).fname;
        fexist = dir(imgname);
        if isempty(fexist)
          unix(sprintf('gunzip -fq %s', [imgname, '.gz']));
        end
      end
    end

    % Find the contrast number for the effect of interest
    numcons = length(spm_original.SPM.xCon);
    count = 1;
    for i = 1:numcons
      if strcmp(spm_original.SPM.xCon(i).name, 'effects of interest')
        break;
      else
        count = count + 1;
      end
    end
    ppi_xY.Ic = count;

    fprintf('contrast %d is the effects of interest\n', count);

    %-Extract BOLD in the seed regions
    disp('---> Extracting Seed BOLDs ...');
    BOLD1 = ppi_regions(ppi_dir, ppi_xY);

    %-Second ROI
    ppi_xY.xyz = roi_center{2}(:);
    ppi_xY.name = roi_name{2};
    ppi_xY.Sess = session_number;
    ppi_xY.def = 'sphere';
    ppi_xY.spec = roi_radius(2);

    roi_info(2).name = ppi_xY.name;
    roi_info(2).center = ppi_xY.xyz;
    roi_info(2).shape = ppi_xY.def;
    roi_info(2).spec = ppi_xY.spec;
    roi_info(2).unit = 'mm';

    %-Extract BOLD in the second ROI
    BOLD2 = ppi_regions(ppi_dir, ppi_xY);

    xY = struct(BOLD1);
    xY(2) = BOLD2;

    disp('---> Generating PhiPI regressors');
    %-Generate PPI regressors
    PPI = effconn_peb_ppi(fullfile(pwd, 'SPM.mat'), 'phipi', xY, ...
      ppi_u, ppi_name, 1, ppi_dir);

    %-Copy batch_stats.mat to ppi folder
    cd(ppi_dir);
    unix(sprintf('cp -a %s .', fullfile(stats_dir, 'batch_stats.mat')));

    %-Save roi information
    save roi_info.mat roi_info;

    %-Update the batch file
    disp('---> Updating batch file');
    load batch_stats.mat;
    tempsess = matlabbatch{1}.spm.stats.fmri_spec.sess(1);
    %-Number of images in a session and its starting number
    num_sess_images = spm_original.SPM.nscan(session_number);
    sess_start_image = start_number(session_number);
    sess_end_image = sess_start_image + num_sess_images -1;
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = pwd;
    matlabbatch{1}.spm.stats.fmri_spec.sess = tempsess;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {};
    if strcmp(fext, '.nii')
      for i = sess_start_image:sess_end_image
        ii = (i+1) - sess_start_image;
        imgfile = [spm_original.SPM.xY.VY(i).fname, ',', num2str(ii)];
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{ii} = imgfile;
      end
    else
      for i = sess_start_image:sess_end_image
        ii = (i+1) - sess_start_image;
        imgfile = [spm_original.SPM.xY.VY(i).fname, ',1'];
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{ii} = imgfile;
      end
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

    %-PPI interaction term
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = ...
      ['PhiPI_', ppi_name];
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = PPI.ppi;

    %-Seed timeseries
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = ...
      ['SeedBOLD1_', roi_name{1}];
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = PPI.Y;

    %-Contrast of task designs
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = ...
      ['SeedBOLD2_', roi_name{2}];
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = PPI.P;

    matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(pwd, 'SPM.mat');

    matlabbatch{3}.spm.stats.con.spmmat{1} = fullfile(pwd, 'SPM.mat');
    matlabbatch{3}.spm.stats.con.consess = {''};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'PhiPI-Interaction';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [1 0 0 0];

    contrastNames{1} = ['PhiPI_', ppi_name];
    contrastVecs{1} = [1 0 0 0];
    numTContrasts = 1;
    save contrasts.mat contrastNames contrastVecs numTContrasts;

    spm('defaults', 'fmri');
    %-Delete old batch file and SPM.mat file
    unix('/bin/rm -rf batch_stats.mat');
    save batch_stats_phipi matlabbatch
    %-Initialize the batch system
    disp('---> Run PPI estimation');
    spm_jobman('initcfg');
    delete(get(0,'Children'));
    %-Run analysis
    spm_jobman('run', 'batch_stats_phipi.mat');

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

cd(currentdir);
disp('------------------------------------------------------------------');
fprintf('Changing back to the directory: %s \n', currentdir);
c     = fix(clock);
disp('==================================================================');
fprintf('PhiPI analysis finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');

diary off;
clear all;
close all;
end
