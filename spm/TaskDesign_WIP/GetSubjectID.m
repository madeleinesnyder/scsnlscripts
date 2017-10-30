function [availablesubj, IDList] = GetSubjectID(subj_list, expr_name)

% subj_list = '../subjectslist_copy.txt';
% expr_name = 'addition_block';

subjects = GetSubjectsList(subj_list);    % Get the list of subjects

n      = 1;
IDList = {};
availablesubj = {};
for cnt = 1:size(subjects,2)

  folder = sprintf('/fs/musk1/20%s/%s/fmri/%s/behavioral/', ...
                      subjects{cnt}(1:2), subjects{cnt}, expr_name);

  files  = dir([folder '/*.edat']);         % Get the edat file name
  
  if(isempty(files))
    disp(sprintf('Edat file not found in folder : %s',folder));
  else
    [a b c d e] = regexpi(files(1).name, '(.*)-\d+\.\w+');
    IDList{n}        = char(e{1});     % Remove the file extension
    availablesubj{n} = subjects{cnt};
    n = n  + 1;
  end
  
end