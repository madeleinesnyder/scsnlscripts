function p = MDS_adata_basic1(C,T,S,pp,L,opts)
cdir = fileparts(mfilename('fullpath'));
p = struct();
if nargin < 2,
    M = 4; 
    J = 2; 
    C = zeros(M,M,J);
    Cmx = 0.95; 
    C(2:end,1,1) = Cmx; 
    C(4,3,2) = Cmx; 
end
if nargin < 2, 
    T = 400;
end
if nargin < 3, 
    S = 1; 
end
M = size(C,1);
J = size(C,3);
  
for sub=1:S,   
    pdef.s(sub).Phi = 1 + rand(pp,L)/100;
    pdef.s(sub).B = 1 + rand(M,pp)/100;
    pdef.s(sub).v = zeros(J,T); 
    
    blk = 10;  
    vv = 1+mod(ceil( (1:T)/blk),J); 
    for j=1:J,
        pdef.s(sub).v(j,:) = vv == j;
    end 
    pdef.s(sub).v(1,:) = 1; 
end
p = setstructfields(pdef,p);
p.C = C;   
p.U = zeros(M,J);

opts.T = T; opts.S = S;
p = MDS_init(p,pp,L,opts);  
for k=1:length(p.s)
   p.s(k).s = rand(size(p.s(k).s)); 
end
p = MDS_fward(p);
for k=1:length(p.s),
    p.s(k).y = p.s(k).yhat + p.s(k).e;
end 
end