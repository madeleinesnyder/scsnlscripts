
function [p2] = MDS_selfcheck(p1)
if nargin < 1, test(); end
p2 = MDS_fward(p1);
cmps(p1,p2);
end
function cmps(p1,p2)
ss = fieldnames(p1);
for k=1:length(ss) 
    s = ss{k};    
    x1 = p1.(s);
    x2 = p2.(s);     
    if isstruct(x1),
        cmps(x1,x2);
    elseif isnumeric(x1) || islogical(x1),
%         disp(s);
        if ~eqc(x1,x2),
            disp('SELF CHECK ERROR IN FIELD:');
            s            
            assert(false);
        end
    end 
end
end
