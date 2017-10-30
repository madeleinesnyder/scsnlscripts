% Compute mean, error, stats for ROI data for each ROI, nevents 
% VMenon 2006-09

% NEW VERSION.  Adds support for signalchange by returning
% allsubjectsdata. -HNguyen 2007-08

%--------------------------------------------------------------------------
function [roimeans, roierrors, roistats, allsubjectsdata] = roi_stats_activation(signals, roi_indices, event_indices)
% default all ROIs and event pairs

%--------------------------------------------------------------------------
% initialize 
roimeans = {}; roierrors = {}; roistats = {}; 
nsessions = size(signals{1}.data_roi_sess_event{1}, 2);

nrois = size(roi_indices, 2); 
if (nrois == 0) 
    nrois = size(signals{1}.roi_name, 2);
    roi_indices = 1:nrois;
end;
nevents = size(event_indices, 2);
if (nevents == 0) 
    nevents = size(signals{1}.event_name{1}, 2); % get events from session 1
    event_indices = 1:nevents;
end;
nsubjects = size(signals, 2); 

%--------------------------------------------------------------------------
% first generate the mean and stderr data 
for roi = 1:nrois
  allsubjectsdata{roi} = []; 
  for subject = 1:nsubjects
    subjectdata = zeros(1, nevents);
      for session = 1:nsessions 
	subjectdata = subjectdata + signals{subject}.data_roi_sess_event{roi}{session}(event_indices);
      end
      allsubjectsdata{roi} = [allsubjectsdata{roi} ; [subjectdata/nsessions]]; % average across sessions
  end
  
  roimeans{roi} = mean(allsubjectsdata{roi}); % average across subjects
  roierrors{roi} = std(allsubjectsdata{roi})/sqrt(nsubjects); % stderr across subjects
  
end

%--------------------------------------------------------------------------
%Calculate ttest and SignRank significance between event pairs for each ROI
for roi = 1:nrois
    for event1 = 1:nevents
        for event2 = 1:nevents
	    [h, p1{event1}{event2}]= ttest(allsubjectsdata{roi}(:,event1), allsubjectsdata{roi}(:,event2));
	    [p2{event1}{event2}, h]= signrank(allsubjectsdata{roi}(:,event1), allsubjectsdata{roi}(:,event2));
	end
    end

  roistats.ttest{roi} = p1;
  roistats.signrank{roi} = p2;
end
