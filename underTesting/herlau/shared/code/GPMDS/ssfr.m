function s1 = ssfr(s1,s2),
ss = fieldnames(s2);



for i=1:length(ss),
    s = ss{i};
%      if s == 'u',
%          342
%      end
     
     
     if isfield(s2,s),
        if ~isfield(s1,s),
            s1.(s) = s2.(s);
        else
            % it is in s2 and in s1. 
            
            if numel(s1.(s)) > 1, 
                if isstruct(s1.(s)(1)), 
                    for k=1:numel(s1.(s)),
                        sc(k) = ssfr(s1.(s)(k),s2.(s)(k));
                        
                    end
                    s1.(s) = sc; 
                else
                    s1.(s) = s2.(s);
                end
            elseif numel(s1.(s)) == 1 && numel(s2.(s)) > 1, 
                %%
                s10 = s1.(s);
                for k=1:numel(s2.(s)),                    
                    s1.(s)(k) = ssfr(s10,s2.(s)(k)); 
                end 
                
            elseif isstruct(s1.(s)),                
                s1.(s) = ssfr(s1.(s),s2.(s));
            else
                s1.(s) = s2.(s);
                
            end
            
            
        end
        
         
     end
     
%      if isfield(s1,s) && isfield(s2, s) && isstruct(s1.(s)),
%         s1.(s) = ssfr(s1.(s),s2.(s));
%     end
    
    
end
% for j=1:length(s1),
% s1
% s2
%     s1 = setstructfields(s1,s2);
% end

end