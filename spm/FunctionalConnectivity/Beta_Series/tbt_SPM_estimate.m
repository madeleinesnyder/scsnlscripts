%-trial by trial SPM estimation
function tbt_SPM_estimate(ConfigFile)

idstr = '$Id: tbt_SPM_estimate.m 2014-08-13 v1$';
warning('off', 'MATLAB:FINITE:obsoleteFunction')
c     = fix(clock);
fname = sprintf('tbt_SPM_estimate-%d_%02d_%02d-%02d_%02d_%02.0f.log',c);
diary(fname);
disp('==================================================================');
fprintf('Trial by trial estimation starts at %d/%02d/%02d %02d:%02d:%02d\n',c);
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

%-Remember the current directory
current_dir = pwd;

subject_list = ReadList(strtrim(paralist.subjectlist_file));
original_spm_data_dir = paralist.preprocessed_dir; %'/fs/musk2';
tbt_spm_data_dir = paralist.tbtdata_dir; %''
pipeline = paralist.pipeline_type;        %'swgcavr';
task_to_split = paralist.taskname_list; %{'Trained_Acc','Untrained_Acc'};
stats_dir = paralist.stats_dir;      %'volume_repair_4cond_swgcavr_whiz_4run_9dur_Feb11';

%==========================================================================
addpath(genpath('/home/fmri/fmrihome/SPM/spm8/toolbox/ArtRepair'));

if ~exist(tbt_spm_data_dir, 'dir')
    mkdir(tbt_spm_data_dir);
end

subjects = ReadList(subject_list);
num_subj = length(subjects);
num_task_in = length(task_to_split);

spm('defaults', 'fmri');
spm_jobman('initcfg');

for isubj = 1:num_subj
    fprintf('~~~~~~~~~Working on subject # %d~~~~~~~~~~~\n', isubj);
    sub_year = ['20' subjects{isubj}(1:2)];

    subj_original_spm_dir = fullfile(original_spm_data_dir, sub_year, subjects{isubj}, 'fmri/stats_spm8/', stats_dir);
    temp_subj_tbt_spm_dir = fullfile(tbt_spm_data_dir, subjects{isubj}, 'fmri/stats_spm8/temp');
    subj_tbt_spm_dir = fullfile(tbt_spm_data_dir, subjects{isubj}, 'fmri/stats_spm8');

    if ~exist(subj_tbt_spm_dir, 'dir')
        mkdir(subj_tbt_spm_dir);
    end

    if ~exist(temp_subj_tbt_spm_dir, 'dir')
        mkdir(temp_subj_tbt_spm_dir);
    end

    load(fullfile(subj_original_spm_dir, 'batch_stats.mat'));

    num_sess = length(matlabbatch{1}.spm.stats.fmri_spec.sess);

    matlabbatch(3) = [];
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = temp_subj_tbt_spm_dir;

    smoothed_data_dir = cell(num_sess, 1);
    for isess = 1:num_sess
        task_design = load(fullfile(subj_original_spm_dir, ['task_design_sess', num2str(isess), '.mat']));
        smoothed_data_dir{isess} = fileparts(matlabbatch{1}.spm.stats.fmri_spec.sess(isess).scans{1});
        unix(sprintf('gunzip -fq %s', fullfile(smoothed_data_dir{isess}, [pipeline 'I.nii.gz'])));
        num_task = length(task_design.names);
        [~, task_in_loc] = ismember(task_to_split, task_design.names);
        task_out_loc = 1:num_task;
        task_out_loc(task_in_loc) = [];
        new_task_names = {};
        new_task_onsets = {};
        new_task_durations = {};
        for itask_in = 1:num_task_in
            num_trial = length(task_design.onsets{task_in_loc(itask_in)});
            for itrial = 1:num_trial
                new_task_names = [new_task_names, {[task_design.names{task_in_loc(itask_in)}, '_', num2str(itrial)]}];
                new_task_onsets = [new_task_onsets, {task_design.onsets{task_in_loc(itask_in)}(itrial)}];
                new_task_durations = [new_task_durations, {task_design.durations{task_in_loc(itask_in)}(1)}];
            end
        end

        sess_name = task_design.sess_name;
        names = [new_task_names, task_design.names(task_out_loc)];
        onsets = [new_task_onsets, task_design.onsets(task_out_loc)];
        durations = [new_task_durations, task_design.durations(task_out_loc)];
        reg_file = task_design.reg_file;
        reg_names = task_design.reg_names;
        reg_vec = task_design.reg_vec;


        new_task_design_file = fullfile(subj_tbt_spm_dir, ['task_design_sess', num2str(isess), '.mat']);
        save(new_task_design_file, 'sess_name', 'names', 'onsets', 'durations', 'reg_file', 'reg_names', 'reg_vec');
        matlabbatch{1}.spm.stats.fmri_spec.sess(isess).multi{1} = new_task_design_file;

    end

    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(temp_subj_tbt_spm_dir, 'SPM.mat');

    save(fullfile(temp_subj_tbt_spm_dir, 'batch_stats.mat'), 'matlabbatch');

    clear matlabbatch;

    delete(get(0,'Children'));

    spm_jobman('run', fullfile(temp_subj_tbt_spm_dir, 'batch_stats.mat'));

    scsnl_art_redo(temp_subj_tbt_spm_dir, pipeline,  subj_tbt_spm_dir, smoothed_data_dir);

    % copy contrasts.mat, task_design, batch_stats
    unix(sprintf('/bin/cp -af %s %s', fullfile(temp_subj_tbt_spm_dir, 'batch_stats*'), subj_tbt_spm_dir));
    % remove temporary stats
    unix(sprintf('/bin/rm -rf %s', temp_subj_tbt_spm_dir));

    for isess = 1:num_sess
        unix(sprintf('gzip -fq %s', fullfile(smoothed_data_dir{isess}, [pipeline 'I.nii.gz'])));
    end

end
