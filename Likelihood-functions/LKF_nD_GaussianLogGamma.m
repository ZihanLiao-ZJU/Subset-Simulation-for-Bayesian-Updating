classdef LKF_nD_GaussianLogGamma
    % Likelihood Function: Gaussian + LogGamma Mixture in N-D
    %
    % Features both skewed and symmetric components
    %
    % Analytical Log-Evidence:
    %   Ndim = 2     → logZ ≈ -8.19
    %   Ndim = 5     → logZ ≈ -20.47
    %   Ndim = 10    → logZ ≈ -40.94
    %   Ndim = 20    → logZ ≈ -81.89
    %   Ndim = 30    → logZ ≈ -122.82
    %
    % Reference:
    %   Beaujean & Caldwell (2013), arXiv:1304.7808

    properties
        Ndim
        Pdis
    end

    methods
        function obj = LKF_nD_GaussianLogGamma(Ndim)
            obj.Ndim = Ndim;
            obj.Pdis = Nataf(repmat(makedist("Uniform", "lower", -30, "upper", 30), Ndim, 1), eye(Ndim));
        end

        function L = EvlLKF(~, theta)
            [ndim, Nsam] = size(theta);
            L = zeros(1, Nsam);

            beta = 1; alpha = 1;
            c_a = 10; c_b = -10;
            mu_c = 10; mu_d = -10;
            sd = 1;

            for j = 1:Nsam
                Log_g_a = alpha * (theta(1,j) - c_a) - exp(theta(1,j) - c_a)/beta - log(gamma(alpha)) - alpha * log(beta);
                Log_g_b = alpha * (theta(1,j) - c_b) - exp(theta(1,j) - c_b)/beta - log(gamma(alpha)) - alpha * log(beta);
                Log_g_ab = max(Log_g_a, Log_g_b);

                Log_n_c = -0.5 * ((theta(2,j) - mu_c)/sd)^2 - log(sqrt(2*pi)) - log(sd);
                Log_n_d = -0.5 * ((theta(2,j) - mu_d)/sd)^2 - log(sqrt(2*pi)) - log(sd);
                Log_n_cd = max(Log_n_c, Log_n_d);

                % Combine for first two dimensions
                logL = log(0.5) + Log_g_ab + log(exp(Log_g_a - Log_g_ab) + exp(Log_g_b - Log_g_ab)) ...
                     + log(0.5) + Log_n_cd + log(exp(Log_n_c - Log_n_cd) + exp(Log_n_d - Log_n_cd));

                % Add Gaussian or LogGamma for other dimensions
                for i = 3:ndim
                    if i <= (ndim + 2)/2
                        logL = logL + (alpha * (theta(i,j) - c_a) - exp(theta(i,j) - c_a)/beta - log(gamma(alpha)) - alpha * log(beta));
                    else
                        logL = logL + (-0.5 * ((theta(i,j) - mu_c)/sd)^2 - log(sqrt(2*pi)) - log(sd));
                    end
                end
                L(j) = logL;
            end
        end
    end
end