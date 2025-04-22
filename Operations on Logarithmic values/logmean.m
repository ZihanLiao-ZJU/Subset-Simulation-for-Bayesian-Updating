function B = logmean(A,dim)
% Logarithmic version of sum
% Syntax:
% -----------------------------------------------------------------------------------
% B = logmean(A,dim)
% -----------------------------------------------------------------------------------
%
% Inputs:
% -----------------------------------------------------------------------------------
% A    : logarithmic matrix
% dim  : returns the mean (in logarithm) of the elements along dimension dim
% -----------------------------------------------------------------------------------
%
% Outputs:
% -----------------------------------------------------------------------------------
% B    : mean (in logarithm) array
% -----------------------------------------------------------------------------------

if strcmpi(dim,'all')
    B = logsum(A,"all")-log(numel(A));
else
    B = logsum(A,dim)-log(size(A,dim));
end
return