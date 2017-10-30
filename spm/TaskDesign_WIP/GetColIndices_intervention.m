function [exprname_col, subj_col,session_col,Block_col,equationAcc_col, equationResp_col, equationCresp_col,equationRT_col, onset_col, expName_col] = GetColIndices(str, expr_name)

inds      = regexpi(str, '\t');

ind           = regexpi(str, 'ExperimentName');
exprname_col  = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'Subject');
subj_col  = length(find(inds < ind)) + 1;

ind           = regexpi(str, 'Session');
session_col  = length(find(inds < ind(1))) + 1;


ind           = regexpi(str, 'Equation.OnsetTime');
onset_col  = length(find(inds < ind(1))) + 1;


ind           = regexpi(str, 'TrialType');
expName_col  = length(find(inds < ind(1))) + 1;



if strcmpi(expr_name, 'addition')
  
  ind           = regexpi(str, 'Equation.Acc');
  equationAcc_col    = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.CRESP');
  equationCresp_col  = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.RESP');
  equationResp_col   = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.RT');
  equationRT_col     = length(find(inds < ind(1))) + 1;
  
  
  ind           = regexpi(str, 'Block');
  Block_col     = length(find(inds < ind)) + 1;


elseif strcmpi(expr_name, 'subtraction')
 
  ind           = regexpi(str, 'Equation.Acc');
  equationAcc_col    = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.CRESP');
  equationCresp_col  = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.RESP');
  equationResp_col   = length(find(inds < ind)) + 1;

  ind           = regexpi(str, 'Equation.RT');
  equationRT_col     = length(find(inds < ind(1))) + 1;
  
  
  ind           = regexpi(str, 'Block');
  Block_col     = length(find(inds < ind)) + 1;


  
end




