function PrintROIResults(matfile)

% % matfile = 'signalchange';
% matfile = 'tscore_average';
% % matfile = 'tscore_percent_voxels';
% % matfile = 'beta_average';

filename = sprintf('roi_%s.tsv', matfile);
fid      = fopen(filename ,'w');

load(matfile);

n = 1;  % Pick the first subject
N = size(signal, 2);                              % Number of subjects
R = length(signal{n}.data_roi_sess_event);        % Number of ROI's
S = size(signal{n}.data_roi_sess_event{R}, 2);    % Number of sessions
E = length(signal{n}.data_roi_sess_event{R}{S});  % Number of events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This loop is to print the header information
fprintf(fid, 'Subject\t');
for r = 1:R   

    roi = signal{n}.roi_name{r};    % ROI name

    for s = 1:S

      sess = sprintf('Sn(%d)',s);   % session number

      for e = 1:E
        event = signal{1}.event_name{s}{e};

        if(strcmpi(matfile, 'signalchange'))
          str = sprintf('%s_%s_%s\t', roi, sess, char(event));
        else
            if(strcmpi(matfile, 'beta_average'))
              [a, b, c, ev] = regexpi(char(event), '(Sn\(\d+\).*)');
              event = char(ev{1});
            end
            str = sprintf('%s_%s\t', roi, char(event)); 
        end
        fprintf(fid, str); 
      end

    end % session loop     

end % ROI loop
fprintf(fid, '\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This loop is to print the actual data
for n = 1:N
  [a, b, c, subj] = regexpi(signal{n}.subject_stats_dir, '(\d+-\d+-\d+\.\d)');  % subject folder
  fprintf(fid, '%s\t',subj{1});
  for r = 1:R   

    for s = 1:S

      for e = 1:E
        event = signal{1}.event_name{s}{e};
        
        fprintf(fid, '%.4f\t', signal{n}.data_roi_sess_event{r}{s}(e)); 
      end
      
    end % session loop     
      
  end % ROI loop
  fprintf(fid, '\n');
end % subject loop

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%