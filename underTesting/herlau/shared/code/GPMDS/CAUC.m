
function auc = CAUC(C,Ct)
C0 = C; Ct0 = Ct; 
%% turn each into vectors. 
C = C0; Ct = Ct0;
C = rd(C); Ct = rd(Ct);
C = C(:); Ct = Ct(:) ~= 0;
C = abs(C);
% C = C + rand(size(C))/1000;
% [Ct, C]
% close all;
% [X,Y,T,auc] = perfcurve(Ct, C, 1);
% subplot(1,3,1);  
% plot(T, X, 'm.-');
% subplot(1,3,2); 
% plot(T, Y, 'r.-');
% subplot(1,3,3); 
% plot(X, Y, 'k.-');
%%
auc = calcAUC(C,Ct);

end
