classdef SeMCS
    % Compatibility options are views of the model's complete G snapshot.
    properties
        IntBay
    end
    properties(Dependent)
        Sampler
        G
        NumSam
        ConPrb
        MaxIte
        Verbose
    end
    properties(Access=private)
        SamplerPrototype
    end
    methods
        function obj=SeMCS(s)
            obj.IntBay=MCMCBase.configure(s.BayStc);obj.Sampler=s;
        end
        function s=get.Sampler(obj)
            s=obj.SamplerPrototype;s.BayStc=obj.IntBay;
        end
        function obj=set.Sampler(obj,s)
            [obj.IntBay,obj.SamplerPrototype]=MCMCBase.adopt_sampler(obj.IntBay,s);
        end
        function v=get.G(obj),v=obj.IntBay.G.run.schedule;end
        function obj=set.G(obj,v),obj.IntBay.G.run.schedule=v;end
        function v=get.NumSam(obj),v=obj.IntBay.G.run.NumSam;end
        function obj=set.NumSam(obj,v),obj.IntBay.G.run.NumSam=v;end
        function v=get.ConPrb(obj),v=obj.IntBay.G.p;end
        function obj=set.ConPrb(obj,v),obj.IntBay.G.p=v;end
        function v=get.MaxIte(obj),v=obj.IntBay.G.run.MaxIte;end
        function obj=set.MaxIte(obj,v),obj.IntBay.G.run.MaxIte=v;end
        function v=get.Verbose(obj),v=obj.IntBay.G.run.Verbose;end
        function obj=set.Verbose(obj,v),obj.IntBay.G.run.Verbose=v;end
        function out=RunIte(obj)
            b=obj.IntBay;
            if ismethod(b,'prepare_run'),b=b.prepare_run(obj.Sampler);end
            % Preserve the original legacy cap and seed-selection rules.
            b.G.run.level_offset=1;
            b.G.run.select_at_cap=true;
            % This is the original Engine construction's configure/adopt step.
            worker=obj;worker.IntBay=MCMCBase.configure(b);worker.Sampler=obj.Sampler;
            out=worker.run_levels();
        end
    end
    methods(Access=private)
        function out=run_levels(obj)
            b=MCMCBase.configure(obj.IntBay);s=obj.Sampler;r=b.G.run;
            validateattributes(r.NumSam,{'numeric'},{'vector','integer','positive','finite'});
            validateattributes(r.MaxIte,{'numeric'},{'scalar','integer','positive','finite'});
            cap=max(1,r.MaxIte-r.level_offset);N=repmat(r.NumSam(1),cap+1,1);
            if numel(r.NumSam)>1,N(:)=0;N(1:numel(r.NumSam))=r.NumSam(:);end
            preset=r.schedule;
            if ~isempty(preset),b.G=merge(b.G,preset{1},true);end
            x=cell(cap,1);y=x;g=cell(cap+1,1);final=x;
            counts=zeros(cap,1);shape=zeros(cap,2);calls=zeros(cap,1);records=cell(cap,1);before=records;after=records;
            converged=false;g{1}=b.G;
            for k=1:cap
                if isfield(b.G,'sampler_record'),b.G=rmfield(b.G,'sampler_record');end
                s.BayStc=b;g{k}=b.G;
                if k==1
                    if ismethod(b,'initial_samples') && (~isfield(r,'initial_mode')||~strcmp(r.initial_mode,'sampler')),[x{k},y{k},calls(k),b]=b.initial_samples(N(k));
                    else,[x{k},y{k},calls(k),b]=s.SamGen(N(k),1);end
                else
                    [x{k},y{k},calls(k),b]=s.SamGen(xs,ys,steps);
                end
                shape(k,:)=[size(y{k},2),size(y{k},3)];counts(k)=prod(shape(k,:));
                final{k}=b.G;
                if isfield(b.G,'sampler_record')
                    records{k}=b.G.sampler_record;before{k}=records{k}.before;after{k}=records{k}.after;
                else
                    snapshot=MCMCBase.profile_snapshot(b,'model_initial_samples',struct());
                    records{k}=struct('kernel','model_initial_samples','before',snapshot,'after',snapshot);
                    before{k}=snapshot;after{k}=snapshot;
                end
                if k+1<=numel(preset)
                    b.G=merge(b.G,preset{k+1},false);extra=0;
                else
                    [y{k},~,b,extra]=b.advance(x{k},y{k},b.G.p);
                end
                if ismethod(b,'after_advance'),[b,post_calls]=b.after_advance(x(1:k),y(1:k),g(1:k));extra=extra+post_calls;end
                calls(k)=calls(k)+extra;g{k+1}=b.G;
                if r.stop_before_seeds && b.stop(),converged=true;break;end
                if k==cap && ~r.select_at_cap,break;end
                if N(k+1)<=0,error('semcs:Schedule','Sample budget exhausted.');end
                [xs,ys,steps]=b.select_seeds(x{k},y{k},N(k+1));
                if r.Verbose,fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n',k,counts(k),shape(k,1),shape(k,2),calls(k));end
                if isempty(xs)||k==cap,break;end
                if b.stop(),converged=true;break;end
            end
            out=struct('x',{x(1:k)},'y',{y(1:k)},'g',{g(1:k+1)}, ...
                'Nsam',counts(1:k),'Nsze',shape(1:k,:),'Ncal',calls(1:k), ...
                'Nite',k,'flg_cvg',converged,'intBay',b,'sampled_final_G',{final(1:k)}, ...
                'engine','SeMCS','schema_version',2,'sampler_records',{records(1:k)}, ...
                'sampler_before',{before(1:k)},'sampler_after',{after(1:k)});
        end
    end
end

function next=merge(old,preset,initial)
if ~isstruct(preset),next=old;next.parameters=preset;return;end
next=overlay(old,preset);
% A proposal schedule does not replace the host's execution policy or itself.
next.run=old.run;
% Initial explicit live state is allowed; later levels retain native carry.
if ~initial
 for key={'lambda','HMC_time'}
  name=key{1};if isfield(old,name),next.(name)=old.(name);end
 end
end
end
function result=overlay(result,patch)
for key=fieldnames(patch)'
 name=key{1};value=patch.(name);
 if isfield(result,name)&&isstruct(result.(name))&&isscalar(result.(name))&&isstruct(value)&&isscalar(value)
  result.(name)=overlay(result.(name),value);
 else,result.(name)=value;end
end
end
