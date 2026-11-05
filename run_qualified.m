function result=run_qualified(caseName,describeOnly)
% Explicit registered profile; describeOnly constructs without sampling.
if nargin<2,describeOnly=false;end
root=fileparts(mfilename('fullpath'));addpath(genpath(root));setup_semcs();
switch char(caseName)
 case 'sus_bayes_eggbox2', target=LKF_2D_Eggbox(); N=3340;
 case 'sus_bayes_loggamma2', target=LKF_nD_GaussianLogGamma(2); N=5975;
 case 'sus_bayes_loggamma5', target=LKF_nD_GaussianLogGamma(5); N=7464;
 case 'sus_bayes_shells5', target=LKF_nD_GaussianShells(5); N=10398;
 otherwise,error('SeMCS:Case','Choose a listed qualified case.');
end
model=SuS_Bay_Nataf(target,'closed_ties_unit_fallback_v1');
sampler=aCS(model); simulation=main(sampler);
simulation.NumSam=N;simulation.MaxIte=10000;simulation.ConPrb=.1;
result=struct('case',char(caseName),'N',N,'profile','closed_ties_unit_fallback_v1','estimator','source_native_evidence','modelG',model.G,'samplerG',sampler.BayStc.G);
if describeOnly,return;end
raw=simulation.RunIte();post=NumSim(raw);logZ=post.Z;
result.logZ=logZ;result.Ncal=sum(raw.Ncal);disp(result.logZ);disp(result.Ncal);
end
