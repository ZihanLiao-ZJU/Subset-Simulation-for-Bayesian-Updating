function p = logGauss(u,varargin)
% Compute log pdf of a Gaussian distribution.
% Syntax:
% -------------------------------------------------------------------------------
% p = logGauss(u)
% p = logGauss(u,mu,sigma)
% -------------------------------------------------------------------------------
% InputS:
% -------------------------------------------------------------------------------
% u     : samples                                                     [Ndim,Nsam]
% mu    : mean vector of Gaussian                                        [Ndim,1]
% sigma : covariance matrix of Gaussian                               [Ndim,Ndim]
% -------------------------------------------------------------------------------
% OutputS:
% -------------------------------------------------------------------------------
% p     : pdf in logrithm scale                                          [1,Nsam]
% -------------------------------------------------------------------------------

if nargin<2
    Ndim = size(u,1);
    mu = zeros(Ndim,1);
    sigma = eye(Ndim);
else
    mu = varargin{1};
    sigma = varargin{2};
    [Ndim,k] = size(mu);
    assert(all(size(sigma)==Ndim) && k==1)   % one mu and one Sigma
end
u = u-mu;
[R,p]= chol(sigma);
if p ~= 0
    error('ERROR: Sigma is not Positive Defined!');
end
Q = R'\u;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = Ndim*log(2*pi)+2*sum(log(diag(R)));   % normalization constant
p = -0.5*(c+q);
end

