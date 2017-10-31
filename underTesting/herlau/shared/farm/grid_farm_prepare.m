function grid_farm_prepare(jobname,pars)
[jdir, simfile,fsep] = jobdir(jobname);
if exist(simfile),
%     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');    
%     disp('Attempted to save to existing file. If you want to do this, please remove directory');
%     disp(jdir);
%     disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');        
%     assert(false);
%     return;     
end 
% atomize pars.
names = cell(1, numel(pars));
% return
for a=1:numel(pars),
%     if a == 37,
%         3
%     end
    pig.parfile = sprintf('%s%sgridpar_%g.mat', jdir, fsep,a);
    pig.name = sprintf('%s%sgridpig_%g.mat', jdir, fsep,a);
    pig.runs = 0; 
    pig.a = a;
    pig.com = struct();
    pig.com.t = tic(); 
    pig.com.message = 'raw initialization';
    pig.com.id = -1;
    pig.submitted = false;
    names{a} = pig.name;
    par = pars{a};
    if isfield(par, 'HPC_ppn')
%         disp('had HPC_ppn');
        pig.ppn = par.HPC_ppn;
    else
%         disp('did not have HPC_ppn');
    end
    
    if isfield(par,'DO_RESET_PIG'), 
        if ~par.DO_RESET_PIG, 
            continue; 
        end
    end 
    save_pig(pig);     
    save(pig.parfile, 'par');
end
save(simfile,'names');    

end