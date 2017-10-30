function tbt_beta_fc(ConfigFile)

idstr = '$Id: tbt_beta_fc.m 2014-08-13 v1$';
warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
fname = sprintf('tbt_beta_fc-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp('==================================================================');
fprintf('Trial-by-Trial Beta Series Task Connectivity starts at %d/%02d/%02d %02d:%02d:%02d\n',c);
fprintf('%s\n', idstr);
disp('==================================================================');
disp(['Current directory is: ',pwd]);
disp('------------------------------------------------------------------');

ConfigFile = strtrim(ConfigFile);
%-Check existence of configuration file
if ~exist(ConfigFile,'file')
  fprintf('Cannot find the configuration file ... \n');
  diary off; return;
end

ConfigFile = ConfigFile(1:end-2);

%-Run the configuration file
eval(ConfigFile);

subject_list = paralist.subjectlist_file;
roi_names = paralist.roiname_list;
task_names = paralist.taskname_list;
mask_file = paralist.mask_file;
tbt_spm_data_dir = paralist.tbtdata_dir;
result_dir = paralist.result_dir;

%==========================================================================
current_dir = pwd;

subjects = ReadList(subject_list);
num_subj = length(subjects);
num_roi = length(roi_names);
num_task = length(task_names);

roi_index = cell(num_roi, 1);

mask_index = spm_read_vols(spm_vol(mask_file));
mask_index = find(mask_index ~= 0);

for iroi = 1:num_roi
  roi_index_tmp = spm_read_vols(spm_vol(fullfile(roi_dir, [roi_names{iroi}, '.nii'])));
  roi_index_tmp = find(roi_index_tmp ~= 0);
  roi_index{iroi} = roi_index_tmp;
end

if ~exist(result_dir, 'dir')
  mkdir(result_dir);
end

for isubj = 1:num_subj

  fprintf('processing subject %s\n\n', subjects{isubj});

  subj_tbt_spm_dir = fullfile(tbt_spm_data_dir, subjects{isubj}, 'fmri/stats_spm8');
  cd(subj_tbt_spm_dir);
  load(fullfile(subj_tbt_spm_dir, 'SPM.mat'));
  col_names = SPM.xX.name;

  subj_mask_index = spm_read_vols(spm_vol(fullfile(subj_tbt_spm_dir, 'mask.img')));
  subj_mask_index = find(subj_mask_index ~= 0);
  subj_mask_index = intersect(subj_mask_index, mask_index);

  for iroi = 1:num_roi
    roi_index{iroi} = intersect(roi_index{iroi}, subj_mask_index);
  end

  num_sess = length(SPM.Sess);

  for itask = 1:num_task
    roi_task_betas = [];
    task_whole_brain_betas = [];
    for isess = 1:num_sess
      iG_exp = ['^Sn\(', num2str(isess), ').', task_names{itask}];
      iG_match = regexpi(col_names, iG_exp);
      iG_match = ~cellfun(@isempty, iG_match);
      if sum(iG_match) == 0
        error('confound columns are not well defined or found');
      else
        iG = find(iG_match == 1);
      end

      task_beta_img_v = SPM.Vbeta(iG);
      task_beta_d = spm_read_vols(task_beta_img_v);
      img_dim = size(task_beta_d);
      task_beta_d = reshape(task_beta_d, prod(img_dim(1:3)), img_dim(4));
      task_beta_d = task_beta_d';

      roi_beta = [];
      for iroi = 1:num_roi
        roi_beta = [roi_beta, mean(task_beta_d(:, roi_index{iroi}), 2)];
      end

      roi_beta = roi_beta - repmat(mean(roi_beta, 1), size(roi_beta, 1), 1);

      roi_task_betas = [roi_task_betas; roi_beta];

      task_beta_d = task_beta_d - repmat(mean(task_beta_d, 1), size(task_beta_d, 1), 1);
      task_whole_brain_betas = [task_whole_brain_betas; task_beta_d];
    end

    roi_task_betas = roi_task_betas - repmat(mean(roi_task_betas, 1), size(roi_task_betas, 1), 1);
    task_whole_brain_betas = task_whole_brain_betas - repmat(mean(task_whole_brain_betas, 1), size(task_whole_brain_betas, 1), 1);

    roi_whole_brain_cov_mtx = roi_task_betas'*task_whole_brain_betas./(size(roi_task_betas, 1)-1);
    roi_var = std(roi_task_betas);
    whole_brain_var = std(task_whole_brain_betas);
    roi_whole_brain_var = roi_var'*whole_brain_var;
    task_corr = roi_whole_brain_cov_mtx./roi_whole_brain_var;
    task_corr_Z = 0.5*log((1+task_corr)./(1-task_corr));

    for iroi = 1:num_roi
      roi_task_corr_Z = reshape(task_corr_Z(iroi, :)', img_dim(1:3));
      task_corr_file = fullfile(result_dir, [roi_names{iroi}, '_', task_names{itask}, '_', subjects{isubj}, '_beta_corr_Z.nii']);
      v = SPM.Vbeta(1);
      v.fname = task_corr_file;
      v.private.dat.fname = task_corr_file;
      spm_write_vol(v, roi_task_corr_Z);
    end
  end
  clear SPM;
end

cd(current_dir);
