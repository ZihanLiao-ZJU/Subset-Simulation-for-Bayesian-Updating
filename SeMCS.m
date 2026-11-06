classdef SeMCS
    % Common sequential sampling, ISD updates and native seed selection.

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
                'InitializeProbability', false, ...
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
            intBay = obj.IntBay;
            sampler = obj.Sampler;
            maxi = obj.MaxIte + options.MaxIteOffset;
            if ~isempty(options.MinimumLevels)
                maxi = max(options.MinimumLevels, maxi);
            end

            if options.InitializeProbability
                intBay.G.p = obj.ConPrb;
            end
            if ismethod(intBay, 'start_run')
                intBay = intBay.start_run(obj.NumSam, obj.MaxIte, obj.ConPrb, obj.G);
            end

            targets = zeros(obj.MaxIte + 1, 1);
            if numel(obj.NumSam) > 1
                targets(1:numel(obj.NumSam)) = obj.NumSam;
            else
                targets(:) = obj.NumSam;
            end
            if ~ismethod(intBay, 'start_run') && ~isempty(obj.G)
                intBay.G = SeMCS.prescribed_state(intBay.G, obj.G{1}, true, ...
                    options.PresetMode, options.PreserveStateFields);
            end
            sampler.BayStc = intBay;

            x = {};
            y = {};
            g = {intBay.G};
            Nsam = [];
            Nsze = zeros(0, 2);
            Ncal = [];
            x_seed = [];
            y_seed = [];
            nc = targets(1);
            steps = 1;
            flg_cvg = false;
            probabilities = [];

            for ite = 1:maxi
                % 1. Generate samples from the current ISD.
                if ismethod(intBay, 'sample_level')
                    [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = ...
                        intBay.sample_level(sampler, obj.NumSam, x_seed, y_seed, steps, ite);
                    Nsam(ite, 1) = numel(y{ite}(4, :));
                    Nsze(ite, :) = [size(x{ite}, 2), 1];
                    g{ite, 1} = intBay.G;
                else
                    if ite == 1
                        [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = sampler.SamGen(nc, steps);
                    else
                        [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = sampler.SamGen(x_seed, y_seed, steps);
                    end
                    Nsam(ite, 1) = nc * steps;
                    Nsze(ite, :) = [nc, steps];
                end

                % 2. Update the ISD. ReSuS constructs a separate candidate domain.
                if ismethod(intBay, 'update')
                    next = intBay.update(y{ite}(4, :), y{ite}(3, :), obj.ConPrb);
                    g{ite + 1, 1} = next.G;
                else
                    if ite + 1 <= numel(obj.G)
                        intBay.G = SeMCS.prescribed_state(intBay.G, obj.G{ite + 1}, false, ...
                            options.PresetMode, options.PreserveStateFields);
                        extra = 0;
                    else
                        [y{ite}, ~, intBay, extra] = intBay.UpdObj(x{ite}, y{ite}, obj.ConPrb);
                    end
                    Ncal(ite) = Ncal(ite) + extra;
                    g{ite + 1, 1} = intBay.G;
                    sampler.BayStc = intBay;
                end

                % The hybrid-domain method checks its candidate before selection.
                if options.StopBeforeSeeds
                    if next.converged()
                        flg_cvg = true;
                        break;
                    elseif numel(x) == maxi
                        break;
                    end
                end

                % 3. Select the next seeds using the original numerical rule.
                if ismethod(intBay, 'select_run_seeds')
                    [x_seed, y_seed, steps, intBay, probability] = ...
                        next.select_run_seeds(x{ite}, y{ite}, intBay, obj.NumSam);
                    probabilities(end + 1) = probability;
                else
                    if targets(ite + 1) <= 0
                        error('SeMCS:ScheduleExhausted', 'Sample budget exhausted before seed selection.');
                    end
                    target = targets(ite + 1);
                    if strcmp(options.SeedFunction, 'SltSed')
                        [x_seed, y_seed] = intBay.SltSed(x{ite}, y{ite}, target);
                    else
                        [x_seed, y_seed] = SedSlt(intBay, x{ite}, y{ite}, target);
                    end
                    nc = size(x_seed, 2);
                    steps = round(target / nc);
                    if obj.Verbose
                        fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n', ...
                            ite, Nsam(ite), Nsze(ite, 1), Nsze(ite, 2), Ncal(ite));
                    end
                end
                sampler.BayStc = intBay;

                % 4. Preserve empty-seed/cap priority over native convergence.
                if ~options.StopBeforeSeeds
                    if nc <= 0 || numel(x) == maxi
                        break;
                    end
                    flg_cvg = intBay.FlgCvg;
                    if flg_cvg
                        break;
                    end
                end
            end

            out = struct('x', {x}, 'y', {y}, 'g', {g}, ...
                'Nsam', Nsam, 'Nsze', Nsze, 'Ncal', Ncal, ...
                'Nite', numel(x), 'flg_cvg', flg_cvg, 'intBay', intBay);
            if strcmp(options.OutputSchema, 'probability_product')
                out.conditional_probabilities = probabilities;
                out.sampler = class(sampler);
                out.profile = intBay.profile;
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
