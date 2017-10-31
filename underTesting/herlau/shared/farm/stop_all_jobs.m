function stop_all_jobs(jobname)
disp('marking all jobs as not submitted');
pigs = get_pigs(jobname);
I = find([pigs.submitted]);
for i=1:length(I),    
    pigs(I(i)).submitted = false;
    save_pig(pigs(I(i)));
end



    


end