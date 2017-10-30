function [exprname_col, subj_col, addacc_col addcresp_col addresp_col addRT_col  ctrlacc_col  ctrlcresp_col  ctrlresp_col  ctrlRT_col  stimulus_col jitter_col onset_col] = GetColIndices(str, expr_name)

inds      = regexpi(str, '\t');

ind           = regexpi(str, 'ExperimentName');
exprname_col  = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'Subject');
subj_col  = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'Stimulus');
stimulus_col  = length(find(inds < ind)) + 1;

jitter_col = [];
onset_col  = [];

if(strcmpi(expr_name, 'addition_block'))
  
  ind           = regexpi(str, 'AddStim.Acc');
  addacc_col    = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'AddStim.CRESP');
  addcresp_col  = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'AddStim.RESP');
  addresp_col   = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'AddStim.RT');
  addRT_col     = length(find(inds < ind)) + 1;

elseif(strcmpi(expr_name, 'subtraction_block'))

  ind           = regexpi(str, 'SubStim.Acc');
  addacc_col    = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'SubStim.CRESP');
  addcresp_col  = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'SubStim.RESP');
  addresp_col   = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'SubStim.RT');
  addRT_col     = length(find(inds < ind)) + 1;

elseif(strcmpi(expr_name, 'addition_event_1') || strcmpi(expr_name, 'addition_event_2'))

  ind           = regexpi(str, 'Equation.Acc');
  addacc_col    = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.CRESP');
  addcresp_col  = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.RESP');
  addresp_col   = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.RT');
  addRT_col     = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'JitterNull');
  jitter_col    = length(find(inds < ind)) + 1;
  
  ind           = regexpi(str, 'Equation.OnsetTime');
  onset_col     = length(find(inds < ind)) + 1;
  
  ctrlacc_col   = [];
  ctrlcresp_col = [];
  ctrlresp_col  = [];
  ctrlRT_col    = [];

  return;
end

ind           = regexpi(str, 'ControlStim.Acc');
ctrlacc_col   = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'ControlStim.CRESP');
ctrlcresp_col = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'ControlStim.RESP');
ctrlresp_col  = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'ControlStim.RT');
ctrlRT_col    = length(find(inds < ind)) + 1;


