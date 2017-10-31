function M = flatmat(M,v,doreshape)
if nargin < 2, v = 0; end
if nargin < 3, doreshape = false; end
if ~iscell(M), 
    J = size(M,3);
    h = ceil(sqrt(J)); 
    w = ceil(J/h);
    ms = cell(h,w); 
    for i=1:J,
        ms{i} = M(:,:,i); 
    end
    M = flatmat(ms,v,doreshape);
    return;
end
sz = getHW(numel(M));
if doreshape && ~all(sz == size(M)),
    M2 = cell(sz); 
    for i=1:numel(M), M2{i} = M{i}; end
    M = flatmat(M2,v,false);
    return;
end
%%
J = cellfun(@isempty, M); J = J(:)';
for i=find(J),
    if isempty(i), break; end
    dm = M{find(~J, 1, 'first')};
    if isempty(dm), dm = 0; end
    M{i} = reshape( 0*mod(1:numel(dm),2)/10,size(dm));
end
  %% 
fn = numel(M); 
sz = size(M{1});
for i=1:numel(M),
    x = M{i}; 
    if size(x,3) > 1, x = flatmat(x, v, true); end 
    if numel(x) < sz, x = zeros(sz); end
    x(:,end+1) = v;
    x(end+1,:) = v;
     
    M{i} = x;     
end
%% check same size. 
s0 = size(M{1});
for j=1:numel(M),
    if any(size(M{j})-s0)~=0,
        M{j} = zeros(s0);
    end
    
end
M = cell2mat(M);

if ~isempty(v), 
    M = M(1:end-1, 1:end-1);
end
end