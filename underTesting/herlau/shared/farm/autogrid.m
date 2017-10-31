function res = autogrid(fname,res,opts,RESUBMIT_CRACHED_JOBS)
dopts.ncores = 1; 
dopts.ON_GRID = true;
dopts.paths = {};
dopts.walltime = 8;
dopts.ppn = 1; % requested processors per node. 
dopts.que = 'compute'; % OR compute.
dopts.maxque = 32; 
dopts.quepause = 1; 
dopts.Irerun = [];
dopts.matlab_version = '910';


dopts.stanford = is_stanford();
if dopts.stanford,
    dopts.que = 'normal'; 
    dopts.username = 'herlau';
    dopts.walltime = 4;
    disp('Running on Stanford cluster');
else
    disp('Running on DTU cluster'); 
end
dopts.quos = []; 
dopts.mem = 4000;


if isfield(opts, 'ncores') && opts.ncores > 1, dopts.quepause = 20; end

if nargin < 3, opts = struct(); end
default_job = ~isfield(opts, 'jobname'); 
if default_job, opts.jobname = fname; end
opts = setstructfields(dopts, opts); 
opts.RES_NAME = res;
opts.paths{end+1} = pwd(); 


if nargin < 4, RESUBMIT_CRACHED_JOBS = false; end
if ~isempty(opts.Irerun), 
    DO_RERUN = true;
    [jdir,simfile,fsep] = jobdir(opts.jobname);
    Irrfile = fullfile(jdir,'Irerun');
    save(Irrfile,'Irerun');    
else
    DO_RERUN = false;
end
    
%%
% fullfile(fname)
% fileparts(fname) 
% mfilename(fname) 
if ~RESUBMIT_CRACHED_JOBS,
    S = functions(str2func(fname));
    fid = fopen(S.file,'r'); s = {};
    while true,
        ds = fgetl(fid);
        if ~ischar(ds), break; 
        else s{end+1} = ds; end
    end
    fclose(fid);
    [s,s_fun] = exparts(s,'grid_fun');
    [s,s_rm] = exparts(s,'grid_rm');

    s_bl = sprintf('pari = 1; pars = cell(size(%s));\n',opts.RES_NAME);
    s_bl = [s_bl, 'S = functions(str2func(mfilename()));\n', ...
        '[f,d] = fileparts(S.file);\n', ...
        'addpath(f);\n',...
        ];
    if DO_RERUN,
        s_bl = [s_bl, sprintf('\n load(''%s'');\n',Irrfile) ];         
    end
    s = insert(s, 'grid_before_loop',sprintf(s_bl));    
    
    s_al = [sprintf('\n grid_farm_prepare(''%s'',pars);\n',opts.jobname), ...
        ''];


    s = insert(s, 'grid_after_loop',sprintf(s_al));
    
    s_wl = ['x = whos();\n' ...
    'pars{pari} = struct(); \n' ...
    'for GLi=1:length(x),\n' ...     
    '    if any(strcmp(x(GLi).name, {''x'',''pars'',''ans''})), continue; end\n', ...
    '    pars{pari}.(x(GLi).name) = eval(x(GLi).name);\n' ...
    'end\n',...
    'pari = pari+1;\n'];
    if DO_RERUN,
        %%
        ss = [opts.RES_NAME,'{'];
        for i=1:length(s_fun),
            sfl = strtrim(s_fun{i});
            if length(sfl) >= length(ss) && strcmp(sfl(1:length(ss)), ss)
                break; 
            end
        end
        rsargs = sfl(length(ss)+1:find(sfl == '}', 1, 'first')-1);
        s_wl = [sprintf('DO_RESET_PIG = Irerun(%s);\n',rsargs), s_wl]; 
    
    end 
    
    s = insert(s, 'grid_within_loop',sprintf(s_wl));
    [jdir,simfile,fsep] = jobdir(opts.jobname);
    if isdir(jdir) && default_job && ~DO_RERUN,
        disp('Warning, REMOVING default job');        
        try
            rmdir(jdir, 's');
        catch err
            disp(err);
        end
    end
    wrfile(collapse(s),[fname, '_grid_prepare.m']);
    rehash  
     
%     return;
    eval(sprintf('%s_grid_prepare',fname));
    jobs = get_farm_jobs(opts.jobname);
    x = load(jobs(1).parfile);
    par= x.par;
    s = structvars(par);
    v = ''; for i=1:size(s,1), v = [v, '\n',s(i,:)]; end
    sfun = collapse(s_fun);
    sfun = [sprintf('function par = %s_grid_fun(par,pig)\n',fname), ...
        v, ...
        '\n', ...
        sfun, ...
        sprintf('par.%s = %s;\n',opts.RES_NAME,opts.RES_NAME),...
        '\n end \n',...
        ''];
    wrfile(sprintf(sfun),[fname, '_grid_fun.m']);
    rehash
else
    % in this case, stop jobs which are submitted and resubmit them.
    disp('Resubmitting halted/crashed jobs.');
    stop_all_jobs(opts.jobname)
end

if ~opts.stanford,
    disp('ssh hpc-fe1');
    disp('showclass compute -r');
end
apply_ensure(fname,sprintf('%s_grid_fun',fname),opts.paths,opts.ON_GRID,opts.ncores,1,opts);
t0 = tic();
pause(5);

MAX_QUERIED_JOBS= opts.maxque;

jobs = get_farm_jobs(opts.jobname);   
Imissing = find([jobs.runs] < 1 & ~[jobs.submitted]);

Imia = zeros(1,length(jobs));
% 
while true,        
    jobs = get_farm_jobs(opts.jobname);     
    
    % check: has all jobs finished running.
    running = sum([jobs.runs] ==  0); 
    if running == 0, 
        x = load(jobs(1).parfile);         
        res = x.par.(opts.RES_NAME);
        for i=1:length(jobs),
            disp(jobs(i).parfile)
            x = load(jobs(i).parfile);
            jj = find(cellfun(@(v)~isempty(v), x.par.(opts.RES_NAME)));            
            res{jj} = x.par.(opts.RES_NAME){jj}; 
            x = [];
            pause(0.1);
            
        end 
        break;
    end
    % check: how many jobs are currently querried?
    [~, nq] = qstat_info(opts);
    if ~isempty(Imissing) && nq < round(MAX_QUERIED_JOBS/ 2),
        % the que is relatively small. begin to feed the que with jobs. 
        v = min([MAX_QUERIED_JOBS-nq,length(Imissing)]);
        I = Imissing(1:v); 
        Imissing = Imissing(v+1:end); 
        % now start the jobs in I.
        fprintf('Queing a new chunk of jobs of length %i; remaining %i\n',v,length(Imissing));
        start_jobs(opts.jobname, I, opts.quepause);        
    end
    
    fprintf('Waiting for job to finish, missing %i jobs, running for %g seconds.\n',running,toc(t0));
    m =[[jobs.runs] ; [jobs.submitted] ; [jobs.a]];
    m(:,m(1,:)==1) = [];
    for i=3:size(m,1),
        for j=1:size(m,2),
            if m(2,j) == 1, s = '!'; else s = ''; end
            if ~any(Imissing == m(i,j)) && m(2,j) ~= 1, continue; end
            fprintf('%i%s ', m(i,j),s);
        end
        fprintf('\n');
    end
    %% find dead jobs.
    if isempty(Imissing),
       jobinf = qstat_info(opts);
        if isempty(jobinf), jobsQR = [];
        else
            jobsQR = [jobinf.job_state];  
        end
        jobsQR = jobinf(jobsQR == 'R' | jobsQR == 'Q');
        JI = arrayfun(@(job)strncmpi(job.Job_Name(4:end), opts.jobname,length(opts.jobname)), jobsQR);
        jobsQR = jobsQR(JI);
        
        missingpigs = setdiff(find([jobs.runs]==0 & [jobs.submitted] == 1), ...
            [jobsQR.pig_id]);
        Imia(missingpigs) = Imia(missingpigs)+1;
        fprintf('Missing: %s\n', sprintf('%g ',missingpigs));   
        
        % those that have been miss
        %ing for two iterations or more.
        deadpigs = find(Imia>1);
        Imia(deadpigs) = 0;  
        if any(deadpigs) && ~any(any([jobsQR.job_state] == 'Q')),
            % mark all these as not submitted.
            for i=1:length(deadpigs),
                jobs(deadpigs(i)).submitted = false;
                save_pig(jobs(deadpigs(i)));            
            end
            %resubmit.  
            Imissing = [Imissing,deadpigs];
            fprintf('MIA: %s\n', sprintf('%g ',deadpigs));   
        end
    end
    
    
    %% print output of last job.
    a = find([jobs.submitted], 1, 'first');    
    if ~isempty(a),
        a = jobs(a).a;
        cout = sprintf('outputfile%i',a);
        s = fullfile(fileparts(jobs(a).parfile), cout);
        ss = sprintf('tail -n 20 %s', s);
        fprintf(sprintf('Result of: tail -n 10 grid/%s/%s\n', fname,cout));
        [v,res] = system(ss);
        fprintf(res);
        fprintf('---\n');
    end
    pause(16);
end
end
function s = insert(s, tag, text)
% replace instances of tag.
I = cellfun(@(ds)~isempty(strfind(ds, tag)), s); 
if ~any(I), return; end
s{I} = text;
end
function wrfile(s,fname),
fid = fopen(fname,'w');
fprintf(fid,'%s',s);
fclose(fid);
end
function text = collapse(s)
text = '';
for i=1:length(s),
    text = sprintf('%s%s\n',text,s{i});    
end
end
function [srest, stag] = exparts(s,tag)
tstart = ['<',tag,'>'];
tend = ['</',tag,'>'];

stag = {};
srest = {};

extracting = false;
for i=1:length(s), 
    if ~extracting,
        if ~isempty(strfind(s{i}, tstart));
            extracting = true; 
        else
            srest{end+1} = s{i};
        end        
    else
        if ~isempty(strfind(s{i}, tend)); 
            extracting = false; 
        else
            di = [find(s{i}=='%',1,'first')-1,length(s{i})];
%             di = length(s{i});
            stag{end+1} = s{i}(1:di(1));
        end
    end
end
end