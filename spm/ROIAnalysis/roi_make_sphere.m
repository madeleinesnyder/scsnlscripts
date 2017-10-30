% This script make only sphere ROIs
% It reads in a configuration file. You can find a template of it as named:
% roi_myrois.m.template
%
% How to run:
% start matlab (ml7spm8) and in matlab command prompt: 
% type: roi_make_sphere('myrois.m')
%
%__________________________________________________________________________
% 2009 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: roi_make_sphere.m 2009-07-15 $
% -------------------------------------------------------------------------

function roi_make_sphere(Config_File)

warning('off', 'MATLAB:FINITE:obsoleteFunction')
disp(['Current directory is: ',pwd]);

% -------------------------------------------------------------------------
% Check existence of the configuration file
% -------------------------------------------------------------------------

if(exist(Config_File,'file')==0)
  fprintf('Cannot find the configuration file ... \n');
  return;
end

addpath /home/fmri/fmrihome/SPM/spm8/toolbox/marsbar

Config_File = strtrim(Config_File);
Config_File = Config_File(1:end-2);
eval(Config_File);
clear Config_File;

if exist(roi_folder,'dir')
  cd(roi_folder);
else
  mkdir(roi_folder);
  cd(roi_folder);
end

disp('Making ROIS ...');

for i=1:length(myroi)
  coords = myroi{i}.coords;
  radius = myroi{i}.radius;
  name = myroi{i}.name;

  roi = maroi_sphere(struct('centre', coords, 'radius', radius));
  roi = label(roi, name);

  n = num2str(i);
  r = num2str(radius);
  x = num2str(coords(1));
  y = num2str(coords(2));
  z = num2str(coords(3));
  
  if length(n) == 1
    n = ['0' n];
  end
  
  filename = [n '-' r 'mm_' name '_' x '_' y '_' z '_roi.mat'];
  fpath    = fullfile(roi_folder,filename); 
  save(fpath, 'roi');
end

disp('Making ROIs is done.');

end