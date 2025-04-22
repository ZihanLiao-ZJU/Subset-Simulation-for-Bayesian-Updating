function B = logcumsum(A,dim)
% Logarithmic version of cumsum
% Syntax:
% -----------------------------------------------------------------------------------
% B = logcumsum(A,dim)
% -----------------------------------------------------------------------------------
%
% Inputs:
% -----------------------------------------------------------------------------------
% A    : logarithmic matrix
% dim  : returns the logarithmic cumulative sum of the elements along dimension dim
% -----------------------------------------------------------------------------------
%
% Outputs:
% -----------------------------------------------------------------------------------
% B    : logarithmic cumulative sum array
% -----------------------------------------------------------------------------------
[Nrow,Ncol] = size(A);
B = zeros(Nrow,Ncol);
if dim==1
    for k=1:Nrow
        B(k,:) = logsum(A(1:k,:),1);
    end
elseif dim==2
    for k=1:Ncol
        B(:,k) = logsum(A(:,1:k),2);
    end
end
return