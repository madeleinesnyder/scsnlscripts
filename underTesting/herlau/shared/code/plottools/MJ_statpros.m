function D = MJ_statpros(D)
J = size(D{1}.C_m{1},3);
tic();
for m=1:length(D),
    for j=1:J,
        %squeeze(mean(cellfun(@(c)std(c( ))), D{k}.C_m),1))
        C_m =D{m}.C_m; M = size(C_m{1},1); 
        for xv = 1:size(C_m,2), %(:,xh,xv)
            for xh = 1:size(C_m,3),
                q = C_m(:,xv,xh);
                for i=1:length(q), 
                    q{i} = q{i}(:,:,j);
                end
                S = length(q); 
                q = cell2mat( reshape(q, [1, 1, S]));
                A1 = []; W1 = [];
                
                OSTAT =[];
                for s=1:S,
                    dq = q(:,:,s); dq = dq(find(~eye(M))); 
                    [~,OSTAT(s,:)] = sort(dq);                    
                end
                %%
                if true,
                    %%
                    %close all;
                    T = 200; 
                    X = rand(T,M*(M-1),S); 
                    for i=1:T
                        for s=1:S
                            X(i,:,s) = randperm(size(X,2)); 
                        end
                    end
                    X = median(X,3)
                    [y,x] = hist(X(:),100)
                    sd = sqrt(var(X(:)));
                    %m = mean(X(:));
                    if false,
                    plot(x,y/trapz(x,y),'k-'); 
                    hold on; xx = linspace(0,M*(M-1)); 
                    plot(xx, normpdf(xx,m,sd),'r-');
                    end
                    
                    mv = median(OSTAT,1);
                    alpha = 0.05; 
                    p1 = norminv(alpha/2,m, sd);
                    p2 = norminv(1-alpha/2,m, sd);
                    
                    %imagesc(OSTAT)
                end
                %%
                for ci=1:size(q,1), 
                    for cj=1:size(q,2), 
                        A1(ci,cj) = ttest(squeeze(q(ci,cj,:)),0,'Alpha',0.05/(M * (M-1)*J^0) );    
                        W1(ci,cj) = sign(mean(q(ci,cj,:))); 
                    end                    
                end
                %OSTAT = OSTAT-mean(OSTAT(:));
                mval = mean(OSTAT(:));
                A2 = []; W2 = []; 
                A3 = []; W3 = []; 
                A4 = []; W4 = []; 
                %%
                
                
                %%
                for i=1:size(OSTAT,2),
                    A2(1,i) = ttest(OSTAT(:,i)-mval,0,'Alpha',0.05/( M * (M-1)*J) ); 
                    W2(1,i) = sign(mean(OSTAT(:,i)-mval));
                    
                    
                    
                    nn = 3; 
                    dx = (OSTAT >= size(OSTAT,2) - nn-1) - (OSTAT <= nn);
                    dx = mean(dx,1);
                    
                    %%
                    W3(1,i) = sign(dx(1,i)) * (abs(dx(1,i)) > 0.2) ;
                    A3(1,i) = W3(1,i) ~= 0; %sign(dx(1,i)) * (abs(dx(1,i)) > 0.2) ;
                    %%
                    A4(1,i) = median(OSTAT(:,i),1) < p1  || median(OSTAT(:,i),1) > p2; 
                    W4(1,i) = median(OSTAT(:,i),1) - mval; 
                    
                    
                end          
                
                v = zeros(M);
                v(~eye(M)) = A2(:); 
                A2 = v; 
                %%
                v = v*0; v(~eye(M)) = W2(:); W2 = v;
               
               
                
                
                
                
                %%
                D{m}.AW1{xv,xh,j}.A = A1;
                D{m}.AW1{xv,xh,j}.W = W1;
                
                D{m}.AW2{xv,xh,j}.A = A2; 
                D{m}.AW2{xv,xh,j}.W = W2; 
                
                D{m}.AW3{xv,xh,j}.A = l2sq(A3); 
                D{m}.AW3{xv,xh,j}.W = l2sq(W3); 
                
                D{m}.AW4{xv,xh,j}.A = l2sq(A4); 
                D{m}.AW4{xv,xh,j}.W = l2sq(W4); 
            end
        end
    end    
end
end
function X = l2sq(A)
m = ceil(sqrt(length(A)));
X = zeros(m);
X(~eye(m)) = A(:);

end