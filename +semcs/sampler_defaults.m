function opts=sampler_defaults(kernel)
opts=struct('SedSav',false,'Lambda',0.6,'AcpOpt',0.44,'NumIte',5,'Adapt',true,'Ntry',20,'InitialTime',pi/4);
if nargin<1,return;end
opts.AcceptArithmetic='log';opts.CoordinateMode='identity';
switch kernel
 case 'aCS',opts.Sigma0Mode='unit';opts.CarryLambda=true;opts.GroupMode='available';opts.AcceptanceRateMode='attempts';opts.AdaptEmptyGroups=false;
 case 'ES',opts.ESMode='coordinate';opts.SeedUpdateMode='each_step';
 case 'CW_aMH',opts.Sigma0Mode='seed_std';opts.PriorMode='gaussian';
end
end
