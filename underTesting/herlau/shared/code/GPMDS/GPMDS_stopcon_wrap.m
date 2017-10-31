function [stats,p,opts] = GPMDS_stopcon_wrap(p,opts,fout,vdex, GPMDS_fun,CONT)

%% this assume directory is set correctly.
pwd()
if ~isdir('../../grid') % save crap in grid directory to ensure not synced. 
    mkdir('../../grid');
end
df = sprintf('_%i',vdex);
fc = fullfile(['../../grid/' fout df '.mat']);
if CONT,
    x = load(fc);
    p = x.p;
    stats = x.stats;
    p.stats = stats;
    opts = ssfr(x.opts,opts);    
end
[stats,p,opts] = GPMDS_fun(p,opts);
save(fc,'p','stats','opts');

end