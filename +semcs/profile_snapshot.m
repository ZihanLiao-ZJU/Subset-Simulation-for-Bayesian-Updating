function out=profile_snapshot(b,kernel,effective)
opts=struct();if isfield(b.G.sampler,kernel),opts=b.G.sampler.(kernel);end
out=struct('kernel',kernel,'method_profile',semcs.getfield(b.G,'method_profile','public_20261102'), ...
 'method_family',semcs.getfield(b.G,'method_family','custom'),'options',opts,'effective',effective);
end
