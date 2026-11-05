classdef main < semcs.LegacyEngine
    methods
        function obj=main(sampler),setup_semcs();obj@semcs.LegacyEngine(sampler);end
        function out=RunIte(obj)
            if ~SuSClosedFallbackPolicy.enabled(obj.IntBay)
                out=RunIte@semcs.LegacyEngine(obj);return
            end
            assert(isa(obj.Sampler,'semcs.aCS'),'SuSClosedFallback:Sampler','This opt-in variant requires the existing common aCS.');
            b=SuSClosedFallbackPolicy.prepare_run(obj.IntBay);
            b.G.run.level_offset=1;b.G.run.select_at_cap=true;
            engine=semcs.Engine(b,obj.Sampler);out=engine.run();
        end
    end
end
