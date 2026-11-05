classdef MCMCBase
    % Fixed settings belong to the sampler; evolving lambda/time belong to G.
    properties
        BayStc
        SedSav = false
        Lambda = 0.6
        AcpOpt = 0.44
        NumIte = 5
        Adapt = true
        Ntry = 20
        InitialTime = pi / 4
        AcceptArithmetic = 'log'
        CoordinateMode = 'identity'
    end

    properties (Dependent)
        Ndim
        Nfun
    end

    methods

        function d = get.Ndim(obj)
            d = obj.BayStc.Ndim;
        end

        function n = get.Nfun(obj)
            n = obj.BayStc.Nfun;
        end

        function obj = MCMCBase(b)
            if nargin > 0
                obj.BayStc = b;
            end
        end

        function obj = apply_profile(obj)
            if ~isempty(obj.BayStc) && ismethod(obj.BayStc, 'sampler_options')
                patch = obj.BayStc.sampler_options(class(obj));

                for name = fieldnames(patch)'
                    obj.(name{1}) = patch.(name{1});
                end
            end
        end

        function opts = configuration_options(obj)
            opts = struct();
            names = properties(obj);

            for j = 1:numel(names)
                name = names{j};
                if any(strcmp(name, {'BayStc', 'Ndim', 'Nfun'}))
                    continue;
                end

                opts.(name) = obj.(name);
            end
        end

        function [b, u, y, n, save, L, Nc, Ns] = prepare(obj, varargin)
            b = obj.BayStc;
            d = b.Ndim;
            L = eye(d);
            if isfield(b.G, 'L')
                L = b.G.L;
            elseif isfield(b.G, 'L_pi')
                L = b.G.L_pi;
            end

            if isscalar(L)
                L = L * eye(d);
            end

            if ~isequal(size(L), [d, d]) || any(~isfinite(L(:))) || rank(L) ~= d
                error('MCMC:Reference', 'Gaussian factor must be finite nonsingular Ndim by Ndim.');
            end

            if numel(varargin) == 2
                Nc = varargin{1};
                Ns = varargin{2};
                if ~isscalar(Nc) || ~isscalar(Ns) || ~isfinite(Nc) || ~isfinite(Ns) || Nc < 1 || Ns < 1 || Nc ~= floor(Nc) || Ns ~= floor(Ns)
                    error('MCMC:Count', 'Counts must be positive integers.');
                end

                if ismethod(b, 'SamGen')
                    u = b.SamGen(Nc);
                else
                    u = L * randn(d, Nc);
                end

                [b, y, n] = obj.evaluate(b, u);
                save = true;
            elseif numel(varargin) == 3
                u = varargin{1};
                y = varargin{2};
                Ns = varargin{3};
                Nc = size(u, 2);
                n = 0;
                save = obj.SedSav;
            else
                error('MCMC:Arguments', 'SamGen(Nchains,Nsteps) or SamGen(U_seed,Y_seed,Nsteps).');
            end

            if Nc < 1 || Ns < 1 || Nc ~= floor(Nc) || Ns ~= floor(Ns) || ~isequal(size(u), [d, Nc]) || ~isequal(size(y), [b.Nfun, Nc])
                error('MCMC:Shape', 'Invalid chain count, step count or packed seed outputs.');
            end
        end

        function [b, y, n] = evaluate(~, b, u)
            b = b.EvlY(u);
            y = [b.EvlLKF; b.EvlPDF; b.Y];
            n = b.Ncal;
            if ~isequal(size(y), [b.Nfun, size(u, 2)])
                error('MCMC:Output', 'Invalid packed target outputs.');
            end
        end

        function r = reject(obj, new, old, p)
            if strcmp(obj.AcceptArithmetic, 'exp')
                r = exp(new - old) < p;
            else
                delta = new - old;
                r = isnan(delta) | delta < log(max(p, realmin));
            end
        end

        function z = to_reference(~, b, u, L, opts)
            if strcmp(opts.CoordinateMode, 'model')
                z = b.to_reference(u);
            else
                z = L \ u;
            end
        end

        function u = from_reference(~, b, z, L, opts)
            if strcmp(opts.CoordinateMode, 'model')
                u = b.from_reference(z);
            else
                u = L * z;
            end
        end

        function value = state(~, b, key, default)
            value = default;
            if isfield(b.G, key) && ~isempty(b.G.(key))
                value = b.G.(key);
            end
        end

    end
end
