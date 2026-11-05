% Adapted from Zihan Liao's public MATLAB MCMC samplers; repository LICENSE applies.
classdef ES < MCMCBase
    % Coordinate compatibility mode or scalar-angle global elliptical slice.
    properties
        ESMode = 'coordinate'
        SeedUpdateMode = 'each_step'
    end

    methods

        function obj = ES(b)
            if nargin < 1
                b = [];
            end

            obj@MCMCBase(b);
            obj = obj.apply_profile();
        end

        function [u, y, Ncal, b] = SamGen(obj, varargin)
            [b, u0, y0, Ncal, save, L, Nc, Ns] = obj.prepare(varargin{:});
            opts = obj.configuration_options();
            coordinate = strcmp(opts.ESMode, 'coordinate');
            coordinate2pi = strcmp(opts.ESMode, 'coordinate_2pi');
            assert(coordinate || coordinate2pi || strcmp(opts.ESMode, 'global'), 'ES:Mode', 'Choose coordinate, coordinate_2pi or global');
            if coordinate || coordinate2pi
                validateattributes(opts.Ntry, {'numeric'}, {'scalar', 'integer', 'positive', 'finite'});
            end

            d = b.Ndim;
            z0 = obj.to_reference(b, u0, L, opts);
            u = zeros(d, Nc, Ns);
            y = zeros(b.Nfun, Nc, Ns);
            v = randn(d, Nc, Ns);
            p = rand(1, Nc, Ns);

            for k = 1:Ns
                if k == 1 && save
                    u(:,:, k) = u0;
                    y(:,:, k) = y0;
                    continue;
                end

                if ~coordinate && ~coordinate2pi
                    assert(all(isfinite(y0(1,:))), 'ES:Seed', 'Global ESS requires finite current log likelihood');
                    angle = 2 * pi * rand(1, Nc);
                    lo = angle - 2 * pi;
                    hi = angle;
                    ids = 1:Nc;
                    z = z0;
                    u(:,:, k) = u0;
                    y(:,:, k) = y0;

                    while ~isempty(ids)
                        proposed = z0(:, ids) .* cos(angle) + v(:, ids, k) .* sin(angle);
                        proposal = obj.from_reference(b, proposed, L, opts);
                        [b, yp, nn] = obj.evaluate(b, proposal);
                        Ncal = Ncal + nn;
                        rej = obj.reject(yp(1,:), y0(1, ids), p(1, ids, k));
                        yes = ~rej;
                        z(:, ids(yes)) = proposed(:, yes);
                        u(:, ids(yes), k) = proposal(:, yes);
                        y(:, ids(yes), k) = yp(:, yes);
                        if ~any(rej)
                            break;
                        end

                        ids = ids(rej);
                        angle = angle(rej);
                        lo = lo(rej);
                        hi = hi(rej);
                        negative = angle < 0;
                        lo(negative) = angle(negative);
                        hi(~negative) = angle(~negative);
                        angle = lo + (hi - lo) .* rand(1, numel(ids));
                    end

                    z0 = z;
                    u0 = u(:,:, k);
                    y0 = y(:,:, k);
                    continue
                end

                if coordinate2pi
                    assert(strcmp(opts.SeedUpdateMode, 'each_step'), 'ES:SeedUpdate', 'coordinate_2pi requires each_step');
                    ids = 1:Nc;
                    lo = zeros(d, Nc);
                    hi = 2 * pi * ones(d, Nc);
                    z = z0;
                    u(:,:, k) = u0;
                    y(:,:, k) = y0;

                    for j = 1:opts.Ntry
                        angle = rand(d, numel(ids)) .* (hi - lo) + lo;
                        proposed = cos(angle) .* z0(:, ids) + sin(angle) .* v(:, ids, k);
                        proposal = obj.from_reference(b, proposed, L, opts);
                        [b, yp, nn] = obj.evaluate(b, proposal);
                        Ncal = Ncal + nn;
                        rej = obj.reject(yp(1,:), y0(1, ids), p(1, ids, k));
                        yes = ~rej;
                        z(:, ids(yes)) = proposed(:, yes);
                        u(:, ids(yes), k) = proposal(:, yes);
                        y(:, ids(yes), k) = yp(:, yes);
                        if ~any(rej)
                            break;
                        end

                        ids = ids(rej);
                        angle = angle(:, rej);
                        lo = lo(:, rej);
                        hi = hi(:, rej);
                        if j == 1
                            lo = angle - 2 * pi;
                            hi = angle;
                        else
                            positive = angle > 0;
                            lo(~positive) = angle(~positive);
                            hi(positive) = angle(positive);
                        end
                    end

                    z0 = z;
                    u0 = u(:,:, k);
                    y0 = y(:,:, k);
                    continue
                end

                R = sqrt(z0 .^ 2 + v(:,:, k) .^ 2);
                ratio = zeros(size(z0));
                nz = R > 0;
                ratio(nz) = z0(nz) ./ R(nz);
                oldangle = acos(max(-1, min(1, ratio)));
                ids = 1:Nc;
                lo = zeros(d, Nc);
                hi = pi * ones(d, Nc);
                z = z0;
                all_accepted = false;
                u(:,:, k) = u0;
                y(:,:, k) = y0;

                for j = 1:opts.Ntry
                    angle = lo + (hi - lo) .* rand(d, numel(ids));
                    proposed = cos(angle) .* R(:, ids);
                    proposal = obj.from_reference(b, proposed, L, opts);
                    [b, yp, nn] = obj.evaluate(b, proposal);
                    Ncal = Ncal + nn;
                    rej = obj.reject(yp(1,:), y0(1, ids), p(1, ids, k));
                    yes = ~rej;
                    z(:, ids(yes)) = proposed(:, yes);
                    u(:, ids(yes), k) = proposal(:, yes);
                    y(:, ids(yes), k) = yp(:, yes);
                    if ~any(rej)
                        all_accepted = true;
                        break;
                    end

                    ids = ids(rej);
                    angle = angle(:, rej);
                    lo = lo(:, rej);
                    hi = hi(:, rej);
                    upper = angle > oldangle(:, ids);
                    lo(~upper) = angle(~upper);
                    hi(upper) = angle(upper);
                end

                if strcmp(opts.SeedUpdateMode, 'each_step') || all_accepted
                    z0 = z;
                    u0 = u(:,:, k);
                    y0 = y(:,:, k);
                end
            end
        end

    end
end
