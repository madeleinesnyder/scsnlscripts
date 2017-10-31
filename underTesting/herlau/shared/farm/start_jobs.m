function start_jobs(jobname, I, quepause)
if nargin < 3, quepause = 1; end

% que the jobs with indices I.
jobname_escaped = strrep(jobname,' ', '\ ');
for i=1:length(I),
    a = I(i);
    com = sprintf('cd grid/%s/ && qsub torquesub_%g.sh', jobname_escaped, a);
    system(com);
    pause(quepause);
%     if cores > 1, pause(20); end
end
end