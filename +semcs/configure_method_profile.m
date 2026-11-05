function b=configure_method_profile(b,family,name)
% Method-owned sampling configuration only: no threshold/budget/stop changes.
if nargin<3,name=[];end
assert(ismethod(b,'sampler_profile'),'semcs:Profile','The model must define sampler_profile(family,name).');
b=semcs.configure(b);spec=b.sampler_profile(family,name);
assert(isstruct(spec)&&isscalar(spec)&&all(isfield(spec,{'family','name','default_sampler','options','evidence'})), ...
 'semcs:Profile','The method provider must return one complete profile spec.');
for field={'family','name','default_sampler'}
 value=spec.(field{1});
 assert((ischar(value)&&isrow(value)&&~isempty(value))||(isstring(value)&&isscalar(value)&&strlength(value)>0), ...
  'semcs:Profile','Profile identifiers must be nonempty text.');
end
assert(isstruct(spec.options)&&isscalar(spec.options)&&isstruct(spec.evidence)&&isscalar(spec.evidence), ...
 'semcs:Profile','Profile options and evidence must be scalar structs.');
for kernel=fieldnames(spec.options)'
 key=kernel{1};patch=spec.options.(key);
 assert(isstruct(patch)&&isscalar(patch),'semcs:Profile','Each kernel options block must be one struct.');
 for option=fieldnames(patch)',b.G.sampler.(key).(option{1})=patch.(option{1});end
end
b.G.method_family=spec.family;b.G.method_profile=spec.name;
b.G.default_sampler=spec.default_sampler;b.G.method_profile_evidence=spec.evidence;
if ismethod(b,'apply_method_profile'),b=b.apply_method_profile(spec);end
end
