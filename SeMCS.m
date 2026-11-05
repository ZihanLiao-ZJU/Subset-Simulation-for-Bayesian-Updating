classdef SeMCS
    % Common sequential sampling loop; each model owns its level policies.
    properties
        Sampler
        IntBay
        NumSam = 1000
        ConPrb = 0.1
        MaxIte = 10000
        G = {}
        Verbose = true
    end
    methods
        function obj=SeMCS(sampler)
            obj.Sampler=sampler;obj.IntBay=sampler.BayStc;
            defaults=obj.IntBay.run_defaults();
            for name=fieldnames(defaults)',obj.(name{1})=defaults.(name{1});end
        end
        function out=RunIte(obj)
            c=obj.IntBay.start_run(obj,obj.Sampler);
            for level=1:c.limit
                [c,r]=c.model.sample_level(c,level);
                c=c.model.advance_level(c,r,level);
                [done,c]=c.model.run_complete(c,r,'before_seeds');
                if done,break;end
                c=c.model.select_run_seeds(c,r,level);
                [done,c]=c.model.run_complete(c,r,'after_seeds');
                if done,break;end
            end
            out=c.model.finish_run(c);
        end
    end
    methods(Static)
        function c=start_context(owner,sampler,limit,mode)
            model=owner.IntBay;targets=zeros(owner.MaxIte+1,1);
            if numel(owner.NumSam)>1,targets(1:numel(owner.NumSam))=owner.NumSam;
            else,targets(:)=owner.NumSam;end
            if ~isempty(owner.G),model.G=SeMCS.prescribed_state(model.G,owner.G{1},true,mode);end
            sampler.BayStc=model;
            c=struct('owner',owner,'model',model,'sampler',sampler,'limit',limit, ...
                'targets',targets,'nc',targets(1),'steps',1,'xs',[],'ys',[], ...
                'x',{{}},'y',{{}},'g',{{model.G}},'Nsam',[],'Nsze',zeros(0,2),'Ncal',[],'complete',false);
        end
        function [c,r]=sample_context(c,level)
            if level==1,[x,y,calls,model]=c.sampler.SamGen(c.nc,c.steps);
            else,[x,y,calls,model]=c.sampler.SamGen(c.xs,c.ys,c.steps);end
            c.model=model;
            r=struct('x',x,'y',y,'calls',calls,'count',c.nc*c.steps,'shape',[c.nc,c.steps]);
            c.x{level,1}=x;c.y{level,1}=y;c.Nsam(level,1)=r.count;
            c.Nsze(level,:)=r.shape;c.Ncal(level,1)=calls;
        end
        function out=pack_output(c)
            out=struct('x',{c.x},'y',{c.y},'g',{c.g},'Nsam',c.Nsam, ...
                'Nsze',c.Nsze,'Ncal',c.Ncal,'Nite',numel(c.x),'flg_cvg',c.complete,'intBay',c.model);
        end
        function next=prescribed_state(old,preset,initial,mode)
            if ~isstruct(preset),next=old;next.parameters=preset;
            elseif strcmp(mode,'replace'),next=preset;
            else,next=SeMCS.overlay(old,preset);end
            if ~initial
                for key={'lambda','HMC_time'}
                    name=key{1};if isfield(old,name),next.(name)=old.(name);end
                end
            end
        end
        function result=overlay(result,patch)
            for key=fieldnames(patch)'
                name=key{1};value=patch.(name);
                if isfield(result,name)&&isstruct(result.(name))&&isstruct(value)
                    result.(name)=SeMCS.overlay(result.(name),value);
                else,result.(name)=value;end
            end
        end
    end
end
