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
            sampler.BayStc = intBay;
            stop_before_seeds = isprop(intBay, 'StopBeforeSeeds') && intBay.StopBeforeSeeds;

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
                if ite == 1
                    [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = sampler.SamGen(nc, steps);
                else
                    [x{ite, 1}, y{ite, 1}, Ncal(ite, 1), intBay] = sampler.SamGen(x_seed, y_seed, steps);
                end
                Nsze(ite, :) = [size(x{ite}, 2), size(x{ite}, 3)];
                Nsam(ite, 1) = prod(Nsze(ite, :));
                g{ite, 1} = intBay.G;

                % 2. Form the next ISD; a prescribed G patches the full state.
                if ite + 1 <= numel(obj.G)
                    intBay.G = SeMCS.prescribed_state(intBay.G, obj.G{ite + 1}, false);
                    extra = 0;
                else
                    [y{ite}, ~, intBay, extra] = intBay.UpdObj(x{ite}, y{ite}, intBay.G.p);
                end
                Ncal(ite) = Ncal(ite) + extra;
                g{ite + 1, 1} = intBay.G;
                sampler.BayStc = intBay;

                % A completed ReSuS domain needs no seeds; ordinary methods
                % retain their native stop after the seed operation.
                if ite == maxi || (stop_before_seeds && intBay.FlgCvg)
                    flg_cvg = intBay.FlgCvg;
                    break;
                end

                % 3. The native method chooses seeds; there is no seed router.
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
