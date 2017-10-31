function p = sample_q(p)
if nargin < 1, MDS(); return; end

M = size(p.C,1);
for k=1:length(p.s),    

    Sc_inv = p.s(k).Sqq_inv + p.s(k).w * p.s(k).w';
    %%
    for i=1:M,
        qq = mvnrnd(zeros(M,1)', Sc_inv); 
        p.s(k).Sw_inv_qq(:,i) = qq;  
    end

end
p.s(k).Sw_inv_qq(:) = -inf;

end