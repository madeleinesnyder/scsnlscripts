% if method == 1 : Use (linear) convolution
% if method == 2 : Use Friston standard convolution.
function y = fMRI_bold_convolve(s,dT, method)
if nargin < 3, method = 1; end
if nargin < 1, 
    s = zeros(1,100); 
    s(20:24) = 0.2; 
%     dt = 1; %    
%method = 2;     
%     if false,
%     s = rand(2,100);
    dT = 0.72; 
%     method = 2;
%     end
    y1 = fMRI_bold_convolve(s,dT, 1);
    y2 = fMRI_bold_convolve(s,dT, 2);
    %
    close all;
    tt = (1:length(s))*dT;
    plot(tt,s,'k-'); hold on;
    plot(tt,y1/trapz(tt,abs(y1))); hold all;
    plot(tt,y2/trapz(tt,abs(y2))); hold all;
    
    legend({'Stimulus $s(t)$', 'SPM 3-base response', 'BOLD convolution'});
    %%     
end

%%
[M,T] = size(s);

if method == 1, 
    %%
    % use standard convolution. 
%     close all;
%     clc
%     s = rand(M,T);
%      Phi_mu = spm_hrf(dT); 
     y = zeros(M,T);
     [Phi] = tue_get_model_hrf_3basis(dT);
     Phi = Phi(1,:);
     L = size(Phi,2);
%      s = s * 0;
%      s(1,50) = 1; 
     
     for m=1:M,
         for t=1:T,
             for l=1:L,
                 if t-l+1 == 0, break; end
                 y(m,t) = y(m,t) + s(m,t-l+1) * Phi(1,l);
             end
         end
     end
     if false,     
         plot(1:T, y(1,:)); hold all;
         plot( 50 - 1 + (1:L), Phi,'.-');

    %      y = conv(s(1,:),Phi(1,end:-1:1), 'same')
    %      plot(y); hold all;
    %      plot(s(1,:))
         i = 60; 
         L = size(Phi,2);
         v = s(1,i-L+1:i);
         v(end:-1:1) * Phi(1,:)'
         y(1,i)
     end
elseif method == 2
    [M,T] = size(s);
    
    if M > 1, 
        y = zeros(M,T);
        for m=1:M,
            y(m,:) = fMRI_bold_convolve(s(m,:),dT, method);
        end
        return;
    end
    addpath(fullfile(fileparts(mfilename('fullpath')), '../DCM_deconvolve'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '../../spm12'));
    345;
    %% use the other method.
%     s = sin( .1*(1:T) );
    dT_max = 0.1;
    SCALE = ceil(dT/dT_max);
    
%     SCALE = 3; 
    %
%     ((1:(T*SCALE) )/SCALE)*dT
    
    dT_fast = dT/SCALE;
    s2 = interp1( (1:T), s, (1:1/SCALE:T) );
%     close all;
%     plot(s, 'o'); hold all;
        
    II = 1:SCALE:size(s2,2);
    
    %%
    
% 2345
    f = @f_fmri_bold_tue;
    go = @go_fmri_bold;
    x0 = zeros(6,1);
%     dt = 0.1; % fix me. 
%     SCALE = 1; 
    L = @Lsq;
    w1 = zeros(1,1);
    w2 = zeros(1,1); % this can't be right.
    b = zeros(1,3);
    T = size(s,2);
    opts = struct();
    opts.method = 2; 
    FD = filter_derivative(f, go,dT_fast,x0,s2,b,w1,w2,opts);    
    [EE,x,yhat] = FD.forward(s2,b,[w1,w2],x0);
       y = yhat(II);
%     plot(yhat)
    
    %     s2 = s2(II)
%     plot(s2, '.')
%     length(s2)


%     FD.forward()
%     345
else
    assert(false);
end


end