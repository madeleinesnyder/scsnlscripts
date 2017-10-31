function [jobs,nqh] = qstat_info(opts)
if opts.stanford,
    [jobs,nqh] = qstat_info_slurm(opts);
else
    [jobs,nqh] = qstat_info_qsub(opts);    
end
34;
if ~isempty(jobs),  
   
    nqh = sum(arrayfun(@(j)strcmp(j.job_state, 'Q'), jobs));
else
    nqh = 0; 
end

end
function [jobs,nqh] = qstat_info_slurm(opts)
% 234
%%
opts.username = 'herlau';
nqh = 0; 
[~,out] = system(sprintf('squeue -u %s --Format=JobId,Partition,State,Name:100',opts.username));
out2 = strtrim(out);
out2 = strsplit(out2, '\n');
out2 = out2(2:end);

jobs = [];
for ii=1:length(out2),
     fields = {'job_state','queue','Job_Name'};  
     v = strsplit(out2{ii});
     
    jobs(ii).job_id = v{1};
%     v 
    state = v{3};
    if any(strcmp(state,{'RUNNING','COMPLETING'})),
        state = 'R';
    elseif any(strcmp(state,{'PENDING','CONFIGURING'})),
        state = 'Q';
    else
        state
        keyboard
    end
    jobs(ii).job_state = state;
    jobs(ii).queue = v{2};
    jobs(ii).Job_Name = v{4};
    pid = jobs(ii).Job_Name;
    jobs(ii).pig_id = str2num(pid(find(pid == '.',1,'last')+1:end));
end
    
end

function [jobs,nqh] =qstat_info_qsub(opts)
[~,out] = system(sprintf('qstat -f %s',opts.que));
out2 = out;
%%
out = out2;  
out = strsplit(out,'\n');
I = [find(strncmpi('Job Id',out,6)),length(out)+1]; 

jobs = [];
for ii=1:length(I)-1,
    
    job = out( I(ii):I(ii+1)-1);
    job = cellfun(@strtrim,job,'UniformOutput',false);
    %%
    fields = {'job_state','queue','Job_Name'};
     
    
    for j=1:length(fields),
        field = fields{j};
        J = find(strncmpi(job,field,length(field)));
        s = job{J};
        jobs(ii).(field) = strtrim(s(find(s == '=',1,'first')+1:end));
        
    end
    jobs(ii).job_id = strtrim(job{1}( find(job{1} == ':',1,'first')+1:find(job{1} == '.',1,'first')-1));
    
    pid = jobs(ii).Job_Name; 
%     pid(find(pid == '.',1,'last')+1:end);    
    jobs(ii).pig_id = str2num(pid(find(pid == '.',1,'last')+1:end));
    
    %%
%     job{:}
    
%      break; 
end
command_qstat = sprintf('qstat %s | egrep "Q|H" | wc',opts.que);
[~,o] = system(command_qstat);

nq = str2num(o);    
if length(nq) < 1, 
    nq = MAX_QUERIED_JOBS;
    disp(command_qstat)
    disp(o)
    [~,qs] = system('qstat');
    disp(qs)
    disp('Error, bad-bad output of qstat');
end
nqh = nq(1); % the current number of querried jobs. 

% find(cellfun(@(s)strncmpi('Job Id',out,, out))

%%
end
