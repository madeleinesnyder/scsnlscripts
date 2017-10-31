function [jdir,simfile,fsep] = jobdir(jobname)
curdir = pwd();
fsep = '/';
if curdir(1) == '/', fsep = '/'; end

griddir = [curdir, fsep,'grid'];
if ~isdir(griddir), mkdir(griddir); end
jdir = [griddir, fsep, jobname];
if ~isdir(jdir), mkdir(jdir); end
simfile =[jdir, fsep,'cursim.mat']; 

jdir(jdir == '\') = '/';
simfile(simfile == '\') = '/';
end