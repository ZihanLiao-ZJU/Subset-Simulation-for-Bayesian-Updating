classdef SeMCS
    % Sequential sampling with native ISD updates and seed selection.
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
            if isprop(obj.IntBay, 'Nsam')
                obj.NumSam = obj.IntBay.Nsam;
            end
            if isfield(obj.IntBay.G, 'p')
                obj.ConPrb = obj.IntBay.G.p;
            end
            if isprop(obj.IntBay, 'MaxIte')
                obj.MaxIte = obj.IntBay.MaxIte;
            end
            if isprop(obj.IntBay, 'Verbose')
                obj.Verbose = obj.IntBay.Verbose;
            end
        end

        function out = RunIte(obj)
            % MaxIte bounds sampled levels; a vector NumSam is a finite plan.
            validateattributes(obj.MaxIte, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
            validateattributes(obj.NumSam, {'numeric'}, {'vector', 'integer', 'positive', 'finite', 'nonempty'});
            validateattributes(obj.ConPrb, {'numeric'}, {'scalar', '>', 0, '<', 1});
            targets = double(obj.NumSam(:));
            maxi = obj.MaxIte;
            if numel(targets) == 1
                targets = repmat(targets, maxi, 1);
            else
                maxi = min(maxi, numel(targets));
            end

            intBay = obj.IntBay;
            sampler = obj.Sampler;
            intBay.G.p = obj.ConPrb;
            if ~isempty(obj.G)
                intBay.G = SeMCS.prescribed_state(intBay.G, obj.G{1}, true);
            end
            if ismethod(intBay, 'start_run')
                intBay = intBay.start_run();
            end
            sampler.BayStc = intBay;
            candidate_domain = ismethod(intBay, 'update');

            x = {};
            y = {};
            g = {};
            Nsam = [];
            Nsze = zeros(0, 2);
            Ncal = [];
            x_seed = [];
            y_seed = [];
            nc = targets(1);
            steps = 1;
            flg_cvg = false;

            for ite = 1:maxi
                % 1. Sample the current ISD and retain the actual chain shape.
                if ismethod(intBay, 'sample_level')
                    [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = ...
                        intBay.sample_level(sampler, targets(ite), x_seed, y_seed, steps, ite);
                elseif ite == 1
                    [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = sampler.SamGen(nc, steps);
                else
                    [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = sampler.SamGen(x_seed, y_seed, steps);
                end
                Nsze(ite, :) = [size(x{ite}, 2), size(x{ite}, 3)];
                Nsam(ite, 1) = prod(Nsze(ite, :));
                g{ite, 1} = intBay.G;

                % 2. Form the next ISD; a prescribed G patches the full state.
                if candidate_domain
                    next = intBay;
                    if ite + 1 <= numel(obj.G)
                        next.G = SeMCS.prescribed_state(next.G, obj.G{ite + 1}, false);
                    else
                        next = intBay.update(y{ite}(4, :), y{ite}(3, :), intBay.G.p);
                    end
                    g{ite + 1, 1} = next.G;
                    if next.converged()
                        flg_cvg = true;
                        break;
                    end
                else
                    if ite + 1 <= numel(obj.G)
                        intBay.G = SeMCS.prescribed_state(intBay.G, obj.G{ite + 1}, false);
                        extra = 0;
                    else
                        [y{ite}, ~, intBay, extra] = intBay.UpdObj(x{ite}, y{ite}, intBay.G.p);
                    end
                    Ncal(ite) = Ncal(ite) + extra;
                    g{ite + 1, 1} = intBay.G;
                    sampler.BayStc = intBay;
                end

                % No next chain is needed at the final allocated level.
                if ite == maxi
                    if ~candidate_domain
                        flg_cvg = intBay.FlgCvg;
                    end
                    break;
                end

                % 3. The native method chooses seeds; there is no seed router.
                if candidate_domain
                    [x_seed, y_seed, steps, intBay] = ...
                        next.select_run_seeds(x{ite}, y{ite}, intBay, targets(ite + 1));
                    g{ite + 1, 1} = intBay.G;
                else
                    [x_seed, y_seed] = SedSlt(intBay, x{ite}, y{ite}, targets(ite + 1));
                    nc = size(x_seed, 2);
                    steps = round(targets(ite + 1) / nc);
                    if obj.Verbose
                        fprintf('Level %d: N=%d, chains=%d, steps=%d, calls=%d\n', ...
                            ite, Nsam(ite), Nsze(ite, 1), Nsze(ite, 2), Ncal(ite));
                    end
                    flg_cvg = intBay.FlgCvg;
                    if flg_cvg || nc <= 0
                        break;
                    end
                end
                sampler.BayStc = intBay;
            end

            % g{1:Nite} are sampled states; g{Nite+1} is the last candidate.
            out = struct('x', {x}, 'y', {y}, 'g', {g}, ...
                'Nsam', Nsam, 'Nsze', Nsze, 'Ncal', Ncal, ...
                'Nite', numel(x), 'flg_cvg', flg_cvg, ...
                'intBay', intBay, 'sampler', class(sampler));
        end
    end

    methods (Static)
        function next = prescribed_state(old, preset, initial)
            if isstruct(preset)
                next = old;
                for key = fieldnames(preset)'
                    name = key{1};
                    next.(name) = preset.(name);
                end
            else
                next = old;
                next.parameters = preset;
            end
            % Scheduled density changes do not reset live proposal adaptation.
            if ~initial
                for key = {'lambda', 'HMC_time'}
                    name = key{1};
                    if isfield(old, name)
                        next.(name) = old.(name);
                    end
                end
            end
        end

    end
end
