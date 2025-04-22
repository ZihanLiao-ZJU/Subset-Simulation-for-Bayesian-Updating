function C = logplus (A,B)
% Logarithmic version of A+B
% Syntax:
% -----------------------------------------------------------------------------------
% C = logplus (A,B)
% -----------------------------------------------------------------------------------
%
% Inputs:
% -----------------------------------------------------------------------------------
% A    : logarithmic matrix
% B    : logarithmic matrix
% -----------------------------------------------------------------------------------
%
% Outputs:
% -----------------------------------------------------------------------------------
% C    : added (in logarithm) array
% -----------------------------------------------------------------------------------

% maximum value of each element
maxAB = max(A,B);
C = maxAB + log(exp(A-maxAB)+exp(B-maxAB));
ind1 = A>0;
ind2 = B>0;
ind3 = isinf(A);
ind4 = isinf(B);
ind5 = ~isreal(A) | ~isreal(B);
ind6 = (ind1 & ind3) | (ind2 & ind4);
ind7 = (~ind1 & ind3) & (~ind2 & ind4);
% A=+inf | B=+inf
C(ind6)=inf;
% A=-inf & B=-inf
C(ind7)=-inf;
% A complex or B complex
C(ind5) = log(exp(A(ind5))+exp(B(ind5)));
end