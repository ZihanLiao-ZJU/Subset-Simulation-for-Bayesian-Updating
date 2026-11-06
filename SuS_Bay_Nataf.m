classdef SuS_Bay_Nataf
    % Subset simulation with Nataf transformation for Bayesian updating

    properties
        % LKF     : Likelihood function object (must implement EvlLKF and Pdis.U2X)
        % Nfun    : Number of outputs per sample (default: 4)
        % Ndim    : Dimension of parameter space
        % X       : Samples in standard normal space [Ndim × Nsam]
        % Y       : Function values [4 × Nsam]
        MethodProfile = 'author_latest'
        LKF
        Nfun = 4;
        Ndim
        Ncal
        X
        Y
    end

    properties
        % G       : Parameters for intermediate likelihood function
        %           [L_i, p_i, L_i1, p_i1, z_i1, z_tmp1, z_sum1]
        G = struct('parameters', [-inf, 0, -inf, 0, -inf, -inf, -inf], 'lambda', 0.6);
    end

    properties(Dependent)
        % FlgCvg  : Convergence flag
        FlgCvg
    end

    methods
        function obj = SuS_Bay_Nataf(LKF, method_profile)
            if nargin < 2 || isempty(method_profile)
                method_profile = 'author_latest';
            end
            % Constructor
            % INPUT:
            %   LKF : likelihood function object
            obj.LKF = LKF;
            obj.Ndim = LKF.Ndim;
            if ~any(strcmp(char(method_profile), {'author_latest', 'paper'}))
                error('method_sampler_profile:Name', 'Choose author_latest or paper.');
            end
            obj.MethodProfile = method_profile;
            obj.G.p = 0.1;
        end

        function obj = EvlY(obj, u)
            % Evaluate function values for input samples
            % INPUT:
            %   u : samples in standard normal space [Ndim × Nsam]
            % OUTPUT:
            %   obj.Y : [likelihood; prior] log-values
            x = U2X(obj, u);
            L = obj.LKF.EvlLKF(x);
            P = logGauss(u);
            obj.X = u;
            obj.Y = [L; P];
            obj.Ncal = size(u, 2);
        end

        function Li = EvlLKF(obj)
            % Evaluate intermediate likelihood function (log indicator)
            y = obj.Y;
            L = y(1,:);
            g = obj.G.parameters;
            L_i = g(3);
            if isinf(L_i)
                Li = zeros(size(L));
            else
                Li = log(double(L > L_i));
            end
        end

        function pi = EvlPDF(obj)
            % Evaluate intermediate prior (log prior)
            y = obj.Y;
            p = y(2,:);
            pi = p;
        end

        function [y,g,obj,N_cal] = UpdObj(obj,~,y,p)
            % Update parameters for intermediate distributions
            % INPUT:
            %   y : function values [4 × Nsam]
            %   p : target conditional probability
            % Prescribed G.p remains the live conditional-probability setting.
            p = obj.G.p;
            g = obj.G.parameters;
            L_i = g(3);
            p_i = g(4);
            z_i = g(5);
            z_tmp = g(6);
            z_sum = g(7);
            L = y(3,:);
            Nsam = size(L,2);
            L_sort = sort(L, "descend");
            Nsed = round(Nsam * p);
            if Nsed > 0
                L_i1 = (L_sort(Nsed) + L_sort(Nsed+1)) / 2;
            else
                L_i1 = L_sort(1);
            end
            p_i1 = p_i + log(sum(L > L_i1) / Nsam);

            z_i1 = p_i + logmean(min(logminus(L,L_i), logminus(L_i1,L_i)), "all");
            z_tmp1 = p_i + logmean(logminus(L,L_i), "all");
            z_sum1 = logplus(logplus(logminus(z_sum,z_tmp), z_i), z_tmp1);

            g = [L_i, p_i, L_i1, p_i1, z_i1, z_tmp1, z_sum1];
            obj.G.parameters = g;
            g = obj.G;
            N_cal = 0;
        end

        function [u_sed, y_sed] = SedSlt(obj, u, y, Nsam)
            % Select seed samples for next iteration
            u = u(1:obj.Ndim,:);
            y = y(1:obj.Nfun,:);
            L = y(3,:);
            g = obj.G.parameters;
            L_i1 = g(3);
            Ns = round(Nsam ./ (1:Nsam));
            Nsed = round(Nsam ./ Ns);
            ind_sed = L > L_i1;
            x_sed_can = u(:, ind_sed);
            y_sed_can = y(:, ind_sed);
            Nsed_can = sum(ind_sed, "all");
            Nsed = Nsed(find(Nsed_can - Nsed >= 0, 1, "last"));
            if Nsed_can >= Nsed
                randind = randperm(Nsed_can, Nsed);
            elseif Nsed_can > 0
                randind = randi(Nsed_can, 1, Nsed);
            else
                randind = [];
            end
            u_sed = x_sed_can(:, randind);
            y_sed = y_sed_can(:, randind);
            if ~isempty(y_sed)
                y_sed = obj.UpdY(y_sed);
            end
        end

        function y = UpdY(obj, y)
            % Update function values by recomputing Li and pi
            y = y(3:end,:);
            obj.Y = y;
            Pi = obj.EvlPDF;
            Li = obj.EvlLKF;
            y = [Li; Pi; y];
        end

        function FlgCvg = get.FlgCvg(obj)
            % Compute convergence flag based on change in thresholds and evidence
            g = obj.G.parameters;
            L_i = g(1);
            L_i1 = g(3);
            z_tmp = g(6);
            z_sum = g(7);
            FlgCvg = exp(logminus(L_i1,L_i) - logplus(L_i1,L_i)) < 1e-1 && ...
                     exp(z_tmp - z_sum) < 1e-1;
        end

        function x = U2X(obj, u)
            % Transform from standard normal to physical space via Nataf
            x = obj.LKF.Pdis.U2X(u);
        end
    end

    methods

        function options = sampler_options(obj, kernel)
            name = obj.MethodProfile;
            if isempty(name)
                name = 'author_latest';
            end
            name = char(name);
            if ~any(strcmp(name, {'author_latest', 'paper'}))
                error('method_sampler_profile:Name', 'Choose author_latest or paper.');
            end

            settings = struct();
            settings.aCS = struct('Sigma0Mode', 'seed_std', 'CarryLambda', false);
            settings.ES = struct('ESMode', 'coordinate', 'Ntry', 20);
            settings.aCS.AcceptArithmetic = 'log';
            settings.aCS.CoordinateMode = 'identity';
            settings.aCS.GroupMode = 'available';
            settings.aCS.AcceptanceRateMode = 'attempts';
            settings.ES.AcceptArithmetic = 'log';
            settings.ES.CoordinateMode = 'identity';
            settings.ES.SeedUpdateMode = 'each_step';
            settings.CW_aMH = struct('AcceptArithmetic', 'log', 'PriorMode', 'gaussian', 'CoordinateMode', 'identity');
            if strcmp(name, 'author_latest')
                settings.aCS.AcceptArithmetic = 'exp';
                settings.ES.AcceptArithmetic = 'exp';
                settings.CW_aMH = struct('AcceptArithmetic', 'exp');
                settings.aCS.GroupMode = 'fixed';
            end

            options = struct();
            if isfield(settings, kernel)
                options = settings.(kernel);
            end
        end

    end
end
