
function matplot(C,Ct,name)
%%
cm = abs(rd(C)); cm = max(cm(:));
imagesc(flatmat(rd(C), 0, true));
axis image; colorbar; 
set(gca, 'clim', [-1,1]*cm);
auc = CAUC(C,Ct);
title(sprintf('%s, AUC=%2.2f', name, auc));
%%
end
