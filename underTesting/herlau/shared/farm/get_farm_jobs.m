function jobs = get_farm_jobs(jobname)
[jdir, simfile] = jobdir(jobname);
load(simfile);
for a=1:length(names)
    try 
        jobs(a) = load(names{a});
    catch error,
       error
       disp(names{a})
       fprintf('could not grab job %g\n', a);
       
%        pigs(a) = []
    end    
end
end 