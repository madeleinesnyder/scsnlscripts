% This is the configuration template file for group analysis
% More information can be found in the HELP section below
% _________________________________________________________________________
% 2009-2010 Stanford Cognitive and Systems Neuroscience Laboratory
%
% $Id: groupstats_config.m.template 2010-01-02 $
% -------------------------------------------------------------------------

% Please specify the individualstats server
paralist.stats_path = '/fs/musk2';

% Please specify the parent folder of the static data
% For YEAR data structure, use the first one
% For NONE YEAR data structure, use the second one
paralist.parent_folder        = [''];
% paralist.parent_folder = 'UCSFBoxer'

% Please specify file name holding subjects to be analyzed
% For two group stats, specify two subject list files. For eg.,
% {'group1.txt', 'group2.txt'}.
paralist.subjlist_file = {'all_post.txt'};  

% Please specify the stats folder path
paralist.stats_folder  = 'stats_spm8_AROMA_nonaggr_comparison_5cond_modeling_choice_and_blank';

% Plese specify the folder path to put analysis results
paralist.output_folder = '/mnt/apricot1_share1/MathFUNDamentals/ICA_AROMA/groupstats/post/ALL_no_mvmnt_criteria';

% Please specify the file holding regressors
% If there is no regressor, comment the first line and uncomment the second
% line
%paralist.reg_file = {'IQ.txt','EQ.txt'};
paralist.reg_file = [''];

% Please specify the folder path holding template batches
paralist.template_path = '/home/fmri/fmrihome/SPM/spm8_scripts/BatchTemplates';


% =========================================================================
% HELP on Configuration Setup
% =========================================================================
%
% A sample configuration file (template) can be found at
% /home/fmri/fmrihome/SPM/spm8_scripts/GroupStats/groupstats_config.m.template
%
% We suggest you rename the configuration file before starting the 
% groupstats.m. The input argument should be your configuration file name.
%
% -------------------------- PARAMETER LIST -------------------------------
%
% paralist.<x>
%
% (01). x <- participant_path: 
% The path where the subject folders of your analysis is present. e.g. this
% could be '/fs/plum1_share1/arithmetic_block/participants/'. 
%                       
% (02). x<- subjlist_file: 
% Name of the text file containing the list of subjects. It is assumed that
% file exists in one of the Matlab search paths. If only one list is
% present, you are using one group analysis. If two lists are present, you
% are using two group analysis. 
% One list example: paralist.subjlist_file = {'list.txt'}
% Two list example: paralist.subjlist_file = {'list1.txt','list2.txt'}
%
% (03). x <- stats_folder: 
% Stats folder name. e.g., 'stats_spm5_arabic'.
%
% (04). x <- output_folder: 
% Folder where the group stats outputs is saved.
%
% (05). x <- reg_file: 
% The .txt file containing the covariate of interest. Could be multiple files
% e.g.  {'regressor1.txt','regressor2.txt', ...}
% 
% (06). x <- template_path:
% The folder path holding template batches. Normally, the path is set
% default. You should NOT change it unless your analysis configuration
% parameters are differet from template. Please use the Matlab batch GUI to
% generate your own batch file and put it in your own folder. The batch
% file name should ALWAYS be 'batch_1group' for one group analysis and
% 'batch_2group' for two group analysis.
%
% =========================================================================
