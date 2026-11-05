% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef aCS < MCMCBase
    % Configurable whitened base scale and reset/carry lambda policy.
    properties
        Sigma0Mode = 'unit'
        CarryLambda = true
        GroupMode = 'available'
        AcceptanceRateMode = 'attempts'
    end

    methods

        function obj = aCS(b)
            if nargin < 1
                b = [];
            end

            obj@MCMCBase(b);
            obj = obj.apply_profile();
        end

        function [u, y, Ncal, b] = SamGen(obj, varargin)
            [b, u0, y0, Ncal, save, L, Nc, Ns] = obj.prepare(varargin{:});
            opts = obj.configuration_options();
            d = b.Ndim;
            u = zeros(d, Nc, Ns);
            y = zeros(b.Nfun, Nc, Ns);
            unit = strcmp(opts.Sigma0Mode, 'unit');
            identity = isequal(L, eye(d));
            mapped = strcmp(opts.CoordinateMode, 'model');
            assert(unit || strcmp(opts.Sigma0Mode, 'seed_std'), 'aCS:Sigma0Mode', 'Choose unit or seed_std');
            if (unit || identity) && ~mapped
                v = reshape(L * reshape(randn(d, Nc, Ns), d, []), d, Nc, Ns);
            else
                v = randn(d, Nc, Ns);
            end

            p = rand(1, Nc, Ns);
            if unit
                sigma0 = 1;
                if mapped
                    z0 = obj.to_reference(b, u0, L, opts);
                else
                    z0 = [];
                end
            else
                if identity && ~mapped
                    z0 = u0;
                else
                    z0 = obj.to_reference(b, u0, L, opts);
                end

                sigma0 = max(std(z0, 0, 2), eps);
            end

            if opts.CarryLambda
                lambda = obj.state(b, 'lambda', opts.Lambda);
            else
                lambda = opts.Lambda;
            end

            if ~isscalar(lambda) || ~isfinite(lambda) || lambda <= 0
                error('aCS:State', 'G.lambda must be finite and positive.');
            end

            b.G.lambda = lambda;
            width = max(ceil(Nc / opts.NumIte), 1);
            if strcmp(opts.GroupMode, 'fixed')
                rounds = opts.NumIte;
            else
                rounds = ceil(Nc / width);
            end

            for group = 1:rounds
                first = min((group - 1) * width + 1, Nc);
                last = min(group * width, Nc);
                ids = first:last;
                if unit
                    sigma = min(lambda, 1);
                    rho = sqrt(1 - sigma ^ 2);
                else
                    sigma = min(lambda * sigma0, 1);
                    rho = sqrt(1 - sigma .^ 2);
                end

                count = 0;
                attempts = 0;

                for k = 1:Ns
                    if k == 1 && save
                        u(:, ids, k) = u0(:, ids);
                        y(:, ids, k) = y0(:, ids);
                        continue;
                    end

                    if unit && ~mapped
                        proposal = rho * u0(:, ids) + sigma * v(:, ids, k);
                    elseif identity && ~mapped
                        proposal = rho .* u0(:, ids) + sigma .* v(:, ids, k);
                    else
                        zp = rho .* z0(:, ids) + sigma .* v(:, ids, k);
                        proposal = obj.from_reference(b, zp, L, opts);
                    end

                    [b, yp, nn] = obj.evaluate(b, proposal);
                    Ncal = Ncal + nn;
                    rej = obj.reject(yp(1,:), y0(1, ids), p(1, ids, k));
                    yes = ~rej;
                    u0(:, ids(yes)) = proposal(:, yes);
                    y0(:, ids(yes)) = yp(:, yes);
                    if mapped || (~unit && ~identity)
                        z0(:, ids(yes)) = zp(:, yes);
                    end

                    u(:, ids, k) = u0(:, ids);
                    y(:, ids, k) = y0(:, ids);
                    count = count + sum(yes);
                    attempts = attempts + numel(ids);
                end

                if opts.Adapt && attempts > 0
                    if strcmp(opts.AcceptanceRateMode, 'chain_steps')
                        rate = count / numel(ids) / Ns;
                    else
                        rate = count / attempts;
                    end

                    lambda = exp(log(lambda) + (rate - opts.AcpOpt) / sqrt(group));
                end
            end

            b.G.lambda = lambda;
        end

    end
end
