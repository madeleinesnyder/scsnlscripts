function search_subjects (participant_path)

currentdir = pwd;

participant_paths = ReadList(participant_path);

npath = length(participant_paths);

dirfilt = '[0-9][0-9]-[0-9][0-9]-[0-9][0-9].[0-9]*';

for i = 1:npath
  fprintf('----->');
  if ~exist(participant_paths{i}, 'dir')
    fprintf('Cannot locate %s \n', participant_paths{i});
    continue;
  end
  fprintf('Searching in: %s \n', participant_paths{i});
  fname_full = [num2str(i), '.txt'];
  unix(sprintf('find %s -maxdepth 1 -name ''%s'' > %s', ...
    participant_paths{i}, dirfilt, fname_full));
  
end

cd(currentdir);

disp('--------------------');
disp('Searching is done');

end