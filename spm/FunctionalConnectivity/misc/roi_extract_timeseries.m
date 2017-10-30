% Extract timeseries for specified ROIs by providing image files, filter
% option and the number of beginning images to be truncated
%__________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: roi_extract_timeseries.m 2010-03-26 $
% -------------------------------------------------------------------------

function [timeseries, roi_name] = ...
  roi_extract_timeseries(roi_folder, roi_list, imgs, filter, NUMTRUNC)

                                
disp('ROI timeseries extraction - Starting');
%-Add path of marsbar functions 
addpath /home/fmri/fmrihome/SPM/spm8/toolbox/marsbar;

%-Assume all images within the same data directory
[datadir, fname, fext] = fileparts(deblank(imgs(1,:)));
checkvar = strfind(fext, '.nii');
if ~isempty(checkvar)
  unix(sprintf('gunzip -fq %s', fullfile(datadir, [fname, '*.gz'])));
else
  unix(sprintf('gunzip -fq %s', fullfile(datadir, '*.img.gz')));
  unix(sprintf('gunzip -fq %s', fullfile(datadir, '*.hdr.gz')));
end
  

%-Truncate selected image files
truncpoint = NUMTRUNC+1;
select_data = imgs(truncpoint:end,:);

%-Get ROIs
if isempty(roi_list)  
  ROI_list = get_roilist(roi_folder);
else
  ROIs = multi_readlist(roi_list);
  numROIs = length(ROIs{1});
  ROI_list = cell(numROIs, 1);
  for i = 1:numROIs
    ROI_list{i} = fullfile(roi_folder, ROIs{1}{i});
  end
end
nrois = length(ROI_list);
if nrois==0
  error('No ROIs specified')
end

%-Get timeseries for each ROI
timeseries = [];
roi_name = cell(nrois,1);
for j = 1:nrois
  roi_obj = maroi(ROI_list{j});
  roi_name{j} = label(roi_obj);
  roi_data_obj = get_marsy(roi_obj, select_data, 'mean');
  roi_ts = summary_data(roi_data_obj);
  % Zero-mean and linear detrend if specified
  if filter 
    roi_ts = roi_ts - mean(roi_ts);
    roi_ts = detrend(roi_ts);
  end
  timeseries = [timeseries; roi_ts'];
end

%-Zip files back
if ~isempty(checkvar)
  unix(sprintf('gzip -fq %s', fullfile(datadir, [fname, '*.nii'])));
else
  unix(sprintf('gzip -fq %s', fullfile(datadir, '*.img')));
  unix(sprintf('gzip -fq %s', fullfile(datadir, '*.hdr')));
end

disp('ROI timeseries extraction - Done');

end