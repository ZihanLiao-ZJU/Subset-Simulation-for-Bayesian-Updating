function b=configure(b)
defaults=struct('NumSam',1000,'MaxIte',10000,'Verbose',true, ...
 'stop_before_seeds',false,'select_at_cap',false,'level_offset',0,'schedule',{{}});
if isprop(b,'Nsam'),defaults.NumSam=b.Nsam;end
if ~isfield(b.G,'p'),b.G.p=0.1;end
for key=fieldnames(defaults)'
 name=key{1};if ~isfield(b.G,'run')||~isfield(b.G.run,name),b.G=semcs.setfield(b.G,['run.' name],defaults.(name));end
end
if ~isfield(b.G,'sampler'),b.G.sampler=struct();end
% Legacy flat writes are consumed once, then removed: one canonical source.
flat=semcs.sampler_defaults();legacy=intersect(fieldnames(flat),fieldnames(b.G.sampler));
for kernel={'aCS','ES','HMC','CW_aMH'}
 name=kernel{1};opts=semcs.sampler_defaults(name);
 if isfield(b.G.sampler,name)
  incoming=b.G.sampler.(name);
  for key=fieldnames(incoming)',opts.(key{1})=incoming.(key{1});end
 end
 for key=legacy',opts.(key{1})=b.G.sampler.(key{1});end
 b.G.sampler.(name)=opts;
end
if ~isempty(legacy),b.G.sampler=rmfield(b.G.sampler,legacy);end
end
