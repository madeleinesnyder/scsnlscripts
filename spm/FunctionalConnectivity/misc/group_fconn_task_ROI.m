% Compute mean, error, t and signrank tests for ROI data 
%_________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: group_fconn_task_ROI.m 2010-03-26 $
% -------------------------------------------------------------------------

function rfconnect_group = group_fconn_task_ROI(rfconnect)

%-Default with all sessions and conditions
%-Use z statistics from Fisher's z transformation of correlations

numsubj = length(rfconnect);
numsess = length(rfconnect{1}.data_sess_event_matrix);
numcond = zeros(numsess,1);
for i = 1:numsess
  numcond(i) = length(rfconnect{1}.data_sess_event_matrix{i});
end

numrois = rfconnect{1}.nrois;

%-Mean, std, ttest and signrank of functional connectivities
rfconnect_group = {};
for sesscnt = 1:numsess
  for condcnt = 1:numcond(sesscnt)
    groupdata = [];
    for subjcnt = 1:numsubj
      tempdata = ...
        rfconnect{subjcnt}.data_sess_event_matrix{sesscnt}{condcnt}{3};
      tempdata = tril(tempdata,-1);
      tempdata = tempdata(tempdata~=0);
      groupdata = [groupdata; tempdata'];
    end
    rfconnect_group{sesscnt}{condcnt}.data = groupdata;
    rfconnect_group{sesscnt}{condcnt}.mean = ...
      tril2full(mean(groupdata), numrois);
    rfconnect_group{sesscnt}{condcnt}.mean_std = ...
      tril2full(std(groupdata)/sqrt(numsubj), numrois);
    [h, p] = ttest(groupdata);
    rfconnect_group{sesscnt}{condcnt}.ttest = tril2full(p, numrois);
    numroipair = size(groupdata, 2);
    signrank_stats = zeros(1,numroipair);
    for i = 1:numroipair
      signrank_stats(i) = signrank(groupdata(:,i));
    end
    rfconnect_group{sesscnt}{condcnt}.signrank = ...
      tril2full(signrank_stats, numrois);
  end
end

end

