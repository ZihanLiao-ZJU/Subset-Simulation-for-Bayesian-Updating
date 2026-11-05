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
        function obj = SeMCS(sampler)
            obj.Sampler = sampler;
            obj.IntBay = sampler.BayStc;

            defaults = obj.IntBay.run_defaults();
            for name = fieldnames(defaults)'
                obj.(name{1}) = defaults.(name{1});
            end
        end

        function out = RunIte(obj)
            % Initialize the model and the run-local sampling history.
            context = obj.IntBay.start_run(obj, obj.Sampler);

            for level = 1:context.limit
                % SamGen: draw initial samples or continue the current chains.
                if ismethod(context.model, 'sample_level')
                    [context, record] = context.model.sample_level(context, level);
                else
                    if level == 1
                        [x, y, calls, model] = context.sampler.SamGen(context.nc, context.steps);
                    else
                        [x, y, calls, model] = context.sampler.SamGen(context.xs, context.ys, context.steps);
                    end
                    context.model = model;

                    record = struct( ...
                        'x', x, ...
                        'y', y, ...
                        'calls', calls, ...
                        'count', context.nc * context.steps, ...
                        'shape', [context.nc, context.steps]);

                    context.x{level, 1} = x;
                    context.y{level, 1} = y;
                    context.Nsam(level, 1) = record.count;
                    context.Nsze(level, :) = record.shape;
                    context.Ncal(level, 1) = calls;
                end

                % UpdObj / update: advance the intermediate distribution.
                context = context.model.advance_level(context, record, level);

                % Each model retains its original seed/stop ordering.
                [done, context] = context.model.run_complete(context, record, 'before_seeds');
                if done
                    break;
                end

                % SltSed / SedSlt / select_seeds: choose seeds for the next level.
                context = context.model.select_run_seeds(context, record, level);
                [done, context] = context.model.run_complete(context, record, 'after_seeds');
                if done
                    break;
                end
            end

            out = context.model.finish_run(context);
        end
    end

    methods(Static)
        function context = start_context(owner, sampler, limit, mode)
            model = owner.IntBay;
            targets = zeros(owner.MaxIte + 1, 1);

            if numel(owner.NumSam) > 1
                targets(1:numel(owner.NumSam)) = owner.NumSam;
            else
                targets(:) = owner.NumSam;
            end

            if ~isempty(owner.G)
                model.G = SeMCS.prescribed_state(model.G, owner.G{1}, true, mode);
            end
            sampler.BayStc = model;

            context = struct( ...
                'owner', owner, ...
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
        end

        function out = pack_output(context)
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
        end

        function next = prescribed_state(old, preset, initial, mode)
            if ~isstruct(preset)
                next = old;
                next.parameters = preset;
            elseif strcmp(mode, 'replace')
                next = preset;
            else
                next = SeMCS.overlay(old, preset);
            end

            % Preserve adapted values when applying a parameter schedule.
            if ~initial
                for key = {'lambda', 'HMC_time'}
                    name = key{1};
                    if isfield(old, name)
                        next.(name) = old.(name);
                    end
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
