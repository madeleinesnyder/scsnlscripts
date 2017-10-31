function [p,ps] = flanker_data(M,J)
cdir = fileparts(mfilename('fullpath'));
x = load(fullfile(cdir,'FlankerData'));
p = x.FlankerData{1,1};
ps = x.FlankerData;

end 