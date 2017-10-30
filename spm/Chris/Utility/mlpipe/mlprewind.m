function pos=mlprewind(pid)
% pos=mlprewind(pid)
%  
% Rewinds a read pipe to the origin

pos=mlpseek(pid,0,-1);