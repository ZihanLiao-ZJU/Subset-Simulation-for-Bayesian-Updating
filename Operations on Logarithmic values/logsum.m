function B = logsum(A,dim)
% Logarithmic version of sum
% Syntax:
% -----------------------------------------------------------------------------------
% B = logsum(A,dim)
% -----------------------------------------------------------------------------------
%
% Inputs:
% -----------------------------------------------------------------------------------
% A    : logarithmic matrix
% dim  : returns the sum (in logarithm) of the elements along dimension dim
% -----------------------------------------------------------------------------------
%
% Outputs:
% -----------------------------------------------------------------------------------
% B    : sum (in logarithm) array
% -----------------------------------------------------------------------------------

ind_rel = A == real(A);
A_sze = size(A);
A_rel = -inf(A_sze);
A_cpx = -inf(A_sze);
A_rel(ind_rel) = A(ind_rel);
A_cpx(~ind_rel) = A(~ind_rel);
% maximum value in dim
if isempty(A_rel)
    maxlogx = max(real(A_cpx),[],dim)+pi*1i;
    minlogx = min(real(A_cpx),[],dim)+pi*1i;
elseif isempty(A_cpx)
    maxlogx = max(A_rel,[],dim);
    minlogx = min(A_rel,[],dim);
else
    maxlogx = max(A_rel,[],dim);
    minlogx = min(real(A_cpx),[],dim)+pi*1i;
end
% calculation
B =  maxlogx + log(sum(exp(A-maxlogx),dim));
% index
ind1 = maxlogx<0;
ind2 = minlogx>0;
ind3 = isinf(maxlogx);
ind4 = isinf(minlogx);
ind5 =  ind1 & ind3;
ind6 =  ind2 & ind4; 
ind7 = ~ind1 & ind3;
% max = -inf
B(ind5) = -inf;
% min = inf
B(ind6) = inf;
% max = inf
B(ind7) = inf;
% eliminate the imagination part
ind_img = imag(B)<=eps;
B(ind_img) = real(B);
return