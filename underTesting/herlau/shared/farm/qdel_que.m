function qdel_que(QUE)
if nargin < 1, 
    QUE = 'k40_interactive'; 
    QUE = 'compute'; 
end
% check if stanford
opts.stanford = is_stanford();
opts.username = 'herlau';
opts.que = 'compute';
234;
%%
if opts.stanford
    %% qdel the normal que.
    jobs = qstat_info(opts);
    % 
    %jobs.job_id
    for j=1:length(jobs),
        id = jobs(j).job_id;
        fprintf('deleting jobs %s\n', id);
        system(sprintf('scancel %s', id));
        pause(0.2); 
    end
    
else
%%
% cd(fullfile(
cd(fileparts(mfilename('fullpath')));
%%
% [s,out] = system('qstat | egrep ''R|Q''');
% out = strsplit(out,'\n'); 
% out = out(1:end-1);

jobs = qstat_info(opts);
I = false(size(jobs));
I = I | cellfun(@(s)strcmp(s,QUE),{jobs.queue});
I = I & ([jobs.job_state] == 'R' | [jobs.job_state] == 'Q');

ids = {jobs(I).job_id}; 
% return

for i=1:length(ids)
%     jname = strsplit(out{i}); 
%     jname = jname{1};
    jname = ids{i};
    fprintf('killing: %s\n',jname);
%     jname = jname(1:find(jname == '.',1,'first')-1);    
    system(sprintf('qdel %s',jname));
end
end

end