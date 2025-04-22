function [x_sed,y_sed] = SedSlt(IntBay,x,y,Nsam)
% select the seeds for next iteration
% Syntax:
% ----------------------------------------------------------------------
% [x_sed,y_sed] = SedSlt(IntL,IntP,x,y,p)
% ----------------------------------------------------------------------
% INPUTS:
% IntL  : Intermediate likelihood function
% IntP  : Intermediate prior distribution
% x     : generated samples from sampler                     [Ndim,Nsam]
% y     : target function value of u                            [6,Nsam]
%        --intermediate likelihood function Li
%        --likelihood function L
%        --target function f
%        --intermediate prior PDF pi
%        --prior PDF p
%        --reference distribution PDF q
% p     : filtering ratio
% ----------------------------------------------------------------------
% OUTPUTS:
% x_sed : seed samples
% y_sed : corresponding function list of x_sed
% ----------------------------------------------------------------------

% else select the seeds with random resampling
% weight for random resampling
lograt = ResRat(IntBay,y);
w_rat = exp(lograt);
Nsed_can = sum(w_rat,"all");
% numbers of required seeds
Ns = round(Nsam./(1:Nsam));
Nsed = round(Nsam./Ns);
Nsed = Nsed(find(Nsed_can-Nsed>=0,1,"last"));
% generate seeds
if Nsed>0
    ind = randsample(size(y,2),Nsed,true,w_rat);
    x_sed = x(:,ind);
    y_sed = y(:,ind);
else
    x_sed = [];
    y_sed = [];
end
if ~isempty(y_sed)
    y_sed = IntBay.UpdY(y_sed);
end
end