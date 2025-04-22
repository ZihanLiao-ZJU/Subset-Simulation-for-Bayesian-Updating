function C = logminus (A,B)
% Logarithmic version of A-B
% Syntax:
% -----------------------------------------------------------------------------------
% C = logminus (A,B)
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
% C    : substracted (in logarithm) array
% -----------------------------------------------------------------------------------
maxAB = max(A,B);
C = maxAB + log(exp(A-maxAB)-exp(B-maxAB));
ind1 = A>0;
ind2 = isinf(A);
ind3 = B>0;
ind4 = isinf(B);
ind5 = ~isreal(A) & isreal(B);
ind6 = ~isreal(B) & isreal(A);
ind7 = ~isreal(A) & ~isreal(B);
% A=+inf
ind8 = ind1 & ind2;
C(ind8)=inf;
% A=-inf & B=-inf
ind9 = ~ind3 & ind4 & ~ind1 & ind2;
C(ind9)=-inf;
% A complex, B real
C(ind5) = pi*1i + logplus(A(ind5)-pi*1i,B(ind5));
% A real, B complex
C(ind6) = logplus(A(ind6),B(ind6)-pi*1i);
% A complex, B complex
C(ind7) = logplus(A(ind7),B(ind7)-pi*1i);
end