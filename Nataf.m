classdef Nataf
    % Nataf transformation class
    %
    % Provides transformations between physical variable space (X) and standard
    % normal space (U), given marginal distributions and correlation structure.
    %
    % Supports both forward (U → X) and inverse (X → U) transformations,
    % with correlation adjustment for Lognormal/Normal combinations.
    %
    % Reference:
    %   Liu & Der Kiureghian (1986), Structural Safety 3(1): 1–12

    properties
        Marg_X     % Marginal distributions (array of distribution objects) [Ndim × 1]
        Rho_X      % Correlation matrix in variable space X                [Ndim × Ndim]
        Rho_U      % Transformed correlation matrix in standard normal U   [Ndim × Ndim]
        RhoL_U     % Cholesky decomposition of Rho_U (lower triangular)    [Ndim × Ndim]
        RhoL_X     % Cholesky decomposition of Rho_X (lower triangular)    [Ndim × Ndim]
        Ndim       % Number of variables
    end

    properties (Constant)
        Nfun = 1;
        LogP = true;  % All outputs are assumed to be in log-domain
    end

    methods
        function obj = Nataf(Marg_X, Rho_X)
            % Constructor
            % INPUT:
            %   Marg_X : array of marginal distributions
            %   Rho_X  : correlation matrix in X space
            obj.Marg_X = Marg_X;
            obj.Rho_X = Rho_X;
            obj.Ndim = size(Marg_X, 1);

            % Validate Rho_X
            [obj.RhoL_X, flg_pos] = chol(Rho_X, 'lower');
            if flg_pos ~= 0
                error('Correlation matrix Rho_X is not positive definite.');
            end

            % Transform to standard normal space
            if all(Rho_X == eye(obj.Ndim))
                obj.Rho_U = eye(obj.Ndim);
            else
                Rho_U = eye(obj.Ndim);
                for i = 1:obj.Ndim
                    for j = i+1:obj.Ndim
                        % Analytical conversions for Normal / Lognormal
                        d1 = Marg_X(i).DistributionName;
                        d2 = Marg_X(j).DistributionName;
                        if strcmpi(d1,'Normal') && strcmpi(d2,'Normal')
                            Rho_U(i,j) = Rho_X(i,j);
                        elseif strcmpi(d1,'Normal') && strcmpi(d2,'Lognormal')
                            cv_j = Marg_X(j).std / Marg_X(j).mean;
                            Rho_U(i,j) = Rho_X(i,j) * cv_j / sqrt(log(1 + cv_j^2));
                        elseif strcmpi(d1,'Lognormal') && strcmpi(d2,'Normal')
                            cv_i = Marg_X(i).std / Marg_X(i).mean;
                            Rho_U(i,j) = Rho_X(i,j) * cv_i / sqrt(log(1 + cv_i^2));
                        elseif strcmpi(d1,'Lognormal') && strcmpi(d2,'Lognormal')
                            cv_i = Marg_X(i).std / Marg_X(i).mean;
                            cv_j = Marg_X(j).std / Marg_X(j).mean;
                            Rho_U(i,j) = log(1 + Rho_X(i,j)*cv_i*cv_j) / ...
                                sqrt(log(1 + cv_i^2) * log(1 + cv_j^2));
                        else
                            % Use 2D Gauss-Legendre integration for nonstandard pairs
                            marg = [Marg_X(i); Marg_X(j)];
                            rho_x = Rho_X(i,j);
                            opt.u_max = [8;8];
                            opt.u_min = [-8;-8];
                            Nint = [200;200];
                            [u,w] = GauIntPot(Nint,opt); u = u'; w = w';
                            fi = (marg(1).icdf(normcdf(u(:,1))) - marg(1).mean) / marg(1).std;
                            fj = (marg(2).icdf(normcdf(u(:,2))) - marg(2).mean) / marg(2).std;
                            coef = fi .* fj .* w;

                            rhoU2X_err = @(rho_u) abs(rho_x - sum( ...
                                coef ./ (2*pi*sqrt(1 - rho_u.^2)) .* ...
                                exp(-1./(2*(1 - rho_u.^2)) .* ...
                                (u(:,1).^2 - 2*rho_u.*u(:,1).*u(:,2) + u(:,2).^2)), 1));
                            Rho_U(i,j) = fminbnd(rhoU2X_err, 0, 1);
                        end
                        Rho_U(j,i) = Rho_U(i,j);
                    end
                end
                Rho_U = min(max(Rho_U, -1), 1);
                obj.Rho_U = Rho_U;
            end

            % Validate Rho_U
            [obj.RhoL_U, flg_pos] = chol(obj.Rho_U, 'lower');
            if flg_pos ~= 0
                error('Transformed correlation matrix Rho_U is not positive definite.');
            end
        end

        function X = U2X(obj, U)
            % Transform from standard normal U-space to physical X-space
            Z = pagemtimes(obj.RhoL_U, U);
            X = zeros(size(U));
            for idim = 1:obj.Ndim
                X(idim,:) = obj.Marg_X(idim).icdf(normcdf(Z(idim,:)));
            end
        end

        function U = X2U(obj, X)
            % Transform from physical X-space to standard normal U-space
            Z = zeros(size(X));
            for idim = 1:obj.Ndim
                Z(idim,:) = norminv(obj.Marg_X(idim).cdf(X(idim,:)));
            end
            U = pagemtimes(inv(obj.RhoL_U), Z);
        end

        function logpdf = EvlPDF(obj, X)
            % Evaluate the log joint PDF of samples in X-space
            ndim = obj.Ndim;
            Nsam = size(X,2);
            Z = zeros(ndim, Nsam);
            logpdf_Z = zeros(ndim, Nsam);
            logpdf_Y = zeros(ndim, Nsam);

            for idim = 1:ndim
                Z(idim,:) = norminv(obj.Marg_X(idim).cdf(X(idim,:)));
                logpdf_Z(idim,:) = logGauss(Z(idim,:));
                logpdf_Y(idim,:) = log(obj.Marg_X(idim).pdf(X(idim,:)));
            end
            logpdf_U = logGauss(Z, zeros(ndim,1), obj.Rho_U);
            logpdf = sum(logpdf_Y,1) - sum(logpdf_Z,1) + logpdf_U;
        end

        function logcdf = EvlCDF(obj, X)
            % Evaluate the log CDF of samples in X-space
            ndim = obj.Ndim;
            Nsam = size(X,2);
            Z = zeros(ndim,Nsam);
            for idim = 1:ndim
                Z(idim,:) = norminv(obj.Marg_X(idim).cdf(X(idim,:)));
            end
            logcdf = logGauss(Z, zeros(ndim,1), obj.Rho_U);
        end

        function X = SamGen(obj, N)
            % Generate N random samples in X-space following the marginal distribution
            U = randn(obj.Ndim, N);
            X = obj.U2X(U);
        end
    end
end
