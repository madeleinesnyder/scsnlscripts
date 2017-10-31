function p = sample_y(p)
for k=1:length(p.s),
    kyI = p.s(k).y_impute == 1;        
    if any(kyI)        
        e = mvnrnd(p.NWe.mu,inv(p.NWe.S_inv),nnz(kyI))';
        p.s(k).e(:,kyI) = e;
        p.s(k).y = p.s(k).yhat + p.s(k).e;
    end
end
end