function GenerateTaskDesign(varargin)
%-------------------------------------------------------------------------------------
% GenerateTaskDesign(subj_list, expr_name, data_file, task_dsgn, optional::rest_exists)
%
% subj_list     : Name of the text file containing the list of subjects.
%                 Specify the full path where the text file is.  If only the file name
%                 is specified, it is assumed that file exists in one of
%                 the matlab search paths.
%  
% expr_name     : Session folder in p21.  For example if your data is in
%                '/fs/plum2_share1/2006/06-05-31.1/fmri/addition_block' 
%                 then 'addition_block' would be the the argument.
%
% data_file     : Text file exported from the e-merge file.
%
% task_dsgn     : Task design file name.
%
% rest_exists   : Optional parameter. 0 if rest does not exists and 1 if rest exists.
%-------------------------------------------------------------------------------------

args      = varargin;
expr_name = varargin{2};

if(strcmpi(expr_name, 'addition_block'))
  GenTaskDesignAddSubBlock(args);

elseif(strcmpi(expr_name, 'subtraction_block'))
  GenTaskDesignAddSubBlock(args);

elseif(strcmpi(expr_name, 'addition_event_1'))
  GenTaskDesignAddSubEvent(args);

elseif(strcmpi(expr_name, 'addition_event_2'))
  GenTaskDesignAddSubEvent(args);

else
  disp('Invalid experiment name');
  return;
end
  

