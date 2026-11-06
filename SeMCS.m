classdef SeMCS
    % Common sequential sampling, native model updates and seed selection.

    properties
        Sampler
        IntBay
        NumSam = 1000
        ConPrb = 0.1
        MaxIte = 10000
        G = {}
        Verbose = true
    end

    properties (Access = private)
        RunOptions
    end

    methods
        function obj = SeMCS(sampler)
            obj.Sampler = sampler;
            obj.IntBay = sampler.BayStc;
            obj.RunOptions = struct( ...
                'MaxIteOffset', -1, ...
                'MinimumLevels', [], ...
                'PresetMode', 'replace', ...
                'PreserveStateFields', {{}}, ...
                'UpdateOutputs', 4, ...
                'SeedFunction', 'SltSed', ...
                'StopBeforeSeeds', false, ...
                'OutputSchema', 'levels');

            if isprop(obj.IntBay, 'RunOptions')
                obj.RunOptions = SeMCS.overlay(obj.RunOptions, obj.IntBay.RunOptions);
            end
            if isprop(obj.IntBay, 'Nsam')
                obj.NumSam = obj.IntBay.Nsam;
            end
            if isfield(obj.IntBay.G, 'p')
                obj.ConPrb = obj.IntBay.G.p;
            end
            for key = {'NumSam', 'ConPrb', 'MaxIte', 'Verbose'}
                name = key{1};
                if isfield(obj.RunOptions, name)
                    obj.(name) = obj.RunOptions.(name);
                end
            end
        end

        function out = RunIte(obj)
            options = obj.RunOptions;
            model = obj.IntBay;
            sampler = obj.Sampler;
            limit = obj.MaxIte + options.MaxIteOffset;
            if ~isempty(options.MinimumLevels)
                limit = max(options.MinimumLevels, limit);
            end

            % The two-output native update uses the evolving G.p.
            if options.UpdateOutputs == 2
                model.G.p = obj.ConPrb;
            end
            if ismethod(model, 'start_run')
                model = model.start_run(obj);
            end

            targets = zeros(obj.MaxIte + 1, 1);
            if numel(obj.NumSam) > 1
                targets(1:numel(obj.NumSam)) = obj.NumSam;
            else
                targets(:) = obj.NumSam;
            end

            if options.UpdateOutputs ~= 1 && ~isempty(obj.G)
                model.G = SeMCS.prescribed_state(model.G, obj.G{1}, true, ...
                    options.PresetMode, options.PreserveStateFields);
            end
            sampler.BayStc = model;

            context = struct( ...
                'owner', obj, ...
                'model', model, ...
                'sampler', sampler, ...
                'limit', limit, ...
                'targets', targets, ...
                'nc', targets(1), ...
                'steps', 1, ...
                'xs', [], ...
                'ys', [], ...
                'x', {{}}, ...
                'y', {{}}, ...
                'g', {{model.G}}, ...
                'Nsam', [], ...
                'Nsze', zeros(0, 2), ...
                'Ncal', [], ...
                'complete', false);
            if strcmp(options.OutputSchema, 'probability_product')
                context.conditional_probabilities = [];
            end

            for level = 1:limit
                % SamGen: draw initial samples or continue the current chains.
                if ismethod(context.model, 'sample_level')
                    [context, record] = context.model.sample_level(context, level);
                    context.g{level, 1} = context.model.G;
                else
                    if level == 1
                        [x, y, calls, context.model] = context.sampler.SamGen(context.nc, context.steps);
                    else
                        [x, y, calls, context.model] = context.sampler.SamGen(context.xs, context.ys, context.steps);
                    end
                    record = struct( ...
                        'x', x, ...
                        'y', y, ...
                        'calls', calls, ...
                        'count', context.nc * context.steps, ...
                        'shape', [context.nc, context.steps]);
                end

                context.x{level, 1} = record.x;
                context.y{level, 1} = record.y;
                context.Nsam(level, 1) = record.count;
                context.Nsze(level, :) = record.shape;
                context.Ncal(level, 1) = record.calls;

                % UpdObj / update: retain each native numerical interface.
                if options.UpdateOutputs == 1
                    context.next = context.model.update(record.values, record.logpdf, obj.ConPrb);
                    context.g{level + 1, 1} = context.next.G;
                else
                    if level + 1 <= numel(obj.G)
                        context.model.G = SeMCS.prescribed_state(context.model.G, obj.G{level + 1}, false, ...
                            options.PresetMode, options.PreserveStateFields);
                        extra = 0;
                    elseif options.UpdateOutputs == 2
                        [~, context.model] = context.model.UpdObj(record.y, context.model.G.p);
                        extra = 0;
                    else
                        [context.y{level}, ~, context.model, extra] = context.model.UpdObj(record.x, record.y, obj.ConPrb);
                    end
                    context.Ncal(level) = context.Ncal(level) + extra;
                    context.g{level + 1, 1} = context.model.G;
                    context.sampler.BayStc = context.model;
                end

                % This native domain update checks convergence before seeds.
                if options.StopBeforeSeeds
                    if context.next.converged()
                        context.complete = true;
                        break;
                    elseif numel(context.x) == limit
                        break;
                    end
                end

                % SltSed / SedSlt / select_seeds: choose next-level seeds.
                if ismethod(context.model, 'select_run_seeds')
                    context = context.model.select_run_seeds(context, record, level);
                else
                    if targets(level + 1) <= 0
                        error('SeMCS:ScheduleExhausted', 'Sample budget exhausted before seed selection.');
                    end
                    target = targets(level + 1);
                    if strcmp(options.SeedFunction, 'SltSed')
                        [context.xs, context.ys] = context.model.SltSed(record.x, context.y{level}, target);
                    else
                        [context.xs, context.ys] = SedSlt(context.model, record.x, context.y{level}, target);
                    end
                    context.nc = size(context.xs, 2);
                    context.steps = round(target / context.nc);
                    if obj.Verbose
                        fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n', ...
                            level, record.count, record.shape(1), record.shape(2), context.Ncal(level));
                    end
                end
                context.sampler.BayStc = context.model;

                if ~options.StopBeforeSeeds
                    if context.nc <= 0 || numel(context.x) == limit
                        break;
                    end
                    context.complete = context.model.FlgCvg;
                    if context.complete
                        break;
                    end
                end
            end

            out = struct( ...
                'x', {context.x}, ...
                'y', {context.y}, ...
                'g', {context.g}, ...
                'Nsam', context.Nsam, ...
                'Nsze', context.Nsze, ...
                'Ncal', context.Ncal, ...
                'Nite', numel(context.x), ...
                'flg_cvg', context.complete, ...
                'intBay', context.model);
            if strcmp(options.OutputSchema, 'probability_product')
                out.conditional_probabilities = context.conditional_probabilities;
                out.sampler = class(context.sampler);
                out.profile = context.model.profile;
            end
        end
    end

    methods (Static)
        function next = prescribed_state(old, preset, initial, mode, retained)
            if nargin < 5
                retained = {};
            end
            if ~isstruct(preset)
                next = old;
                next.parameters = preset;
            elseif strcmp(mode, 'replace')
                next = preset;
            else
                next = SeMCS.overlay(old, preset);
            end

            % Keep the actually adapted proposal values in later schedules.
            if ~initial
                for key = {'lambda', 'HMC_time'}
                    name = key{1};
                    if isfield(old, name)
                        next.(name) = old.(name);
                    end
                end
            end

            for key = retained
                name = key{1};
                if ~isfield(next, name)
                    next.(name) = old.(name);
                end
            end
        end

        function result = overlay(result, patch)
            for key = fieldnames(patch)'
                name = key{1};
                value = patch.(name);
                if isfield(result, name) && isstruct(result.(name)) && isstruct(value)
                    result.(name) = SeMCS.overlay(result.(name), value);
                else
                    result.(name) = value;
                end
            end
        end
    end
end
