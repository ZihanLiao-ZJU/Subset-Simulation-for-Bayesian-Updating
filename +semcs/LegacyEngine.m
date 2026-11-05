classdef LegacyEngine
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
        function obj=LegacyEngine(s)
            obj.IntBay=semcs.configure(s.BayStc);obj.Sampler=s;
        end
        function s=get.Sampler(obj)
            s=obj.SamplerPrototype;s.BayStc=obj.IntBay;
        end
        function obj=set.Sampler(obj,s)
            [obj.IntBay,obj.SamplerPrototype]=semcs.adopt_sampler(obj.IntBay,s);
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
            % The old RunIte loop samples at most MaxIte-1 levels.
            b.G.run.level_offset=1;
            b.G.run.select_at_cap=true;
            engine=semcs.Engine(b,obj.Sampler);out=engine.run();
        end
    end
end
