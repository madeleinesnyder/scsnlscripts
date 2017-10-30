% Calculate roi functional connectivities within subjects
%__________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: subject_fconn_task_ROI.m 2010-03-26 $
% -------------------------------------------------------------------------

function sfconnect = subject_fconn_task_ROI(roi_folder, roi_list, ...
  subject_stats_dir, FILT, NUMTRUNC, TR, is_similar_multisession, global_ts)

sfconnect = {};

sfconnect.subject_stats_dir = subject_stats_dir;

%-Load design information
spm_file = fullfile(subject_stats_dir, 'SPM.mat');
if ~exist(spm_file, 'file')
  fprintf('Cannot find SPM.mat in %s \n', subject_stats_dir);
  return;
end
clear SPM;
load(spm_file);
nsess = length(SPM.Sess);
nconds = length(SPM.Sess(1).U);
sfconnect.nsess = nsess;
sfconnect.nconds = nconds;

%-Get images for each session
nscanvec = [0 SPM.nscan];
for numsession = 2:(nsess+1)
  session = numsession-1;
  startimg = sum(nscanvec(1:session))+1;
  endimg = sum(nscanvec(1:numsession));
  
  %-Extract ROI timeseries, one roi per row
  %-No truncation of timeseries right now, truncate it later
  [TS, roi_name] = roi_extract_timeseries(roi_folder, roi_list, ...
    SPM.xY.P(startimg:endimg, :), 0, 0);
  
  %-Adjusting global time series
  if global_ts == 1
    disp('Regressing out global time series ...');
    XReg = [ones(size(TS,2),1), SPM.xGX.rg(startimg:endimg)];
    TS = ((eye(size(TS,2)) - XReg*pinv(XReg'*XReg)*XReg')*TS')';
  end
      
  if FILT == 1
    disp('Detrending ...');
    TS = (detrend(TS'))';
  end
  
  sfconnect.roi_name = roi_name;
  nrois = length(roi_name);
  sfconnect.nrois = nrois;
  
  all_task_timepts = [];
  for cond = 1:nconds
    
    e_s = (session-1)*(nconds+1)+cond;
    if (is_similar_multisession)
      s = session;
      e = cond;
    else
      s = 1;
      e = e_s;
    end
    
    sfconnect.TS{s} = TS;
    
    % select appropriate timepoints for jth condition
    onsets      = SPM.Sess(session).U(cond).ons;
    durations   = SPM.Sess(session).U(cond).dur;
    onset_times = sec2scans(onsets, TR)+1;
    onset_durs  = sec2scans(durations, TR);
    
    %-Do the truncation here
    task_timepts = [];
    for k = 1:length(onset_times)
      task_timepts = [task_timepts, ...
                        onset_times(k)+NUMTRUNC:onset_times(k)+onset_durs(k)-1];
      all_task_timepts = [all_task_timepts, onset_times(k):onset_times(k)+onset_durs(k)-1];
      
    end
    
    TS_task = TS(:,task_timepts)';
    numtp = size(TS_task,1);
    
    %add by hieu to save the task timepoints, so we can extract
    %only the task timeseries
    sfconnect.task_timepts{s}{e} = task_timepts;
    sfconnect.TS_task{s}{e} = TS_task;

    % functional connectivity for jth cond in ith sess
    [rho_matrix, p_matrix] = corr(TS_task);
    z_matrix = rhoToZTransform(rho_matrix, numtp);
    sfconnect.data_sess_event_matrix{s}{e} = {rho_matrix, p_matrix, z_matrix};
    for roi1 = 1:nrois
      for roi2 = 1:nrois
        sfconnect.rhodata_rois_session_event{roi1}{roi2}{s}(e) = ...
          rho_matrix(roi1, roi2);
        sfconnect.pdata_rois_session_event{roi1}{roi2}{s}(e) = ...
          p_matrix(roi1, roi2);
        sfconnect.zdata_rois_session_event{roi1}{roi2}{s}(e) = ...
          z_matrix(roi1, roi2);
      end
    end
    sfconnect.event_name{s}{e} = SPM.Sess(session).U(cond).name{1};
    
  end
  

  %-All conditions
  sfconnect.task_timepts{s}{e+1} = [sfconnect.task_timepts{s}{1:nconds}];
  TS_task = TS'; % all points
  numtp = size(TS_task,1);
  sfconnect.TS_task{s}{e+1} = TS_task;
  [rho_matrix, p_matrix] = corr(TS_task);
  z_matrix = rhoToZTransform(rho_matrix, numtp);
  sfconnect.data_sess_event_matrix{s}{e+1} = {rho_matrix, p_matrix, z_matrix};
  for roi1 = 1:nrois
    for roi2 = 1:nrois
      sfconnect.rhodata_rois_session_event{roi1}{roi2}{s}(e+1) = ...
        rho_matrix(roi1, roi2);
      sfconnect.pdata_rois_session_event{roi1}{roi2}{s}(e+1) = ...
        p_matrix(roi1, roi2);
      sfconnect.zdata_rois_session_event{roi1}{roi2}{s}(e+1) = ...
        z_matrix(roi1, roi2);
    end
  end
  sfconnect.event_name{s}{e+1} = 'all_conditions'; 
  
  %-Rest
  sfconnect.event_name{s}{e+2} = 'rest';
  rest_timepts = setdiff(1:length(TS), all_task_timepts);
  
  if(~isempty(rest_timepts))
    TS_rest      = TS(:,rest_timepts)';
    numtp = size(TS_rest,1);
    [rho_matrix,p_matrix] = corr(TS_rest);
    z_matrix              = rhoToZTransform(rho_matrix, numtp);
  else
    TS_rest    = [];
    rho_matrix = zeros(nrois);
    p_matrix   = zeros(nrois);
    z_matrix   = zeros(nrois);
  end

  sfconnect.task_timepts{s}{e+2} = rest_timepts;
  sfconnect.TS_task{s}{e+2}      = TS_rest;
  sfconnect.data_sess_event_matrix{s}{e+2} = {rho_matrix, p_matrix, z_matrix};
  for roi1 = 1:nrois
    for roi2 = 1:nrois
      sfconnect.rhodata_rois_session_event{roi1}{roi2}{s}(e+2) = ...
        rho_matrix(roi1, roi2);
      sfconnect.pdata_rois_session_event{roi1}{roi2}{s}(e+2) = ...
        p_matrix(roi1, roi2);
      sfconnect.zdata_rois_session_event{roi1}{roi2}{s}(e+2) = ...
        z_matrix(roi1, roi2);
    end
  end  
end
end

function scan = sec2scans(sec, TR)

% converts secs to scans with 'floor'
% note: 0th sec = 0th scan, so add 1 when dealing with indices

scan = floor(sec./TR);

end

function z_matrix = rhoToZTransform(rho_matrix, numtp)

z_matrix = 0.5*log((1 + rho_matrix)./(1 - rho_matrix))*sqrt(numtp-3);

end