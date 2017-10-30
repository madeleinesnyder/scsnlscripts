function position=mlptell(pid)
% position=mlptell(pid)
%
% Returns current position indicator of pipe pid 

position=mlpmex(pid,'tell');

