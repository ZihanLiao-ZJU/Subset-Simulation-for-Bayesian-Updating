function [b,prototype]=adopt_sampler(b,s)
% Import explicitly supplied kernel options, keeping model/run/live G state.
assert(~isa(s,'handle'),'semcs:HandleSamplerUnsupported', ...
    'Use a value-class sampler; a handle sampler cannot be detached without mutating its caller.');
b=semcs.configure(b);
if isa(s,'semcs.MCMCBase')
    opts=s.configuration_options();
    kernel=semcs.kernel_name(s);
    for key=fieldnames(opts)'
        b.G.sampler.(kernel).(key{1})=opts.(key{1});
    end
    b.G.sampler.actual_kernel=kernel;
    prototype=s.unbound_prototype();
else
    % Value-class custom samplers retain the writable BayStc/SamGen protocol.
    prototype=s;prototype.BayStc=[];
end
end
