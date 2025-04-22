classdef LKF_nD_GaussianShells
    % Likelihood Function: Two concentric Gaussian shells in N-D space
    %
    % Used to evaluate posterior sampling methods on multimodal problems
    %
    % Analytical Log-Evidence:
    %   Ndim = 2     → logZ ≈ -1.75
    %   Ndim = 5     → logZ ≈ -5.67
    %   Ndim = 10    → logZ ≈ -14.59
    %   Ndim = 20    → logZ ≈ -36.09
    %   Ndim = 30    → logZ ≈ -60.13
    %
    % Reference:
    %   Feroz & Hobson (2008), MNRAS, 384(2):449–463

    properties
        Ndim
        Pdis
    end

    methods
        function obj = LKF_nD_GaussianShells(Ndim)
            obj.Ndim = Ndim;
            obj.Pdis = Nataf(repmat(makedist("Uniform", "lower", -6, "upper", 6), Ndim, 1), eye(Ndim));
        end

        function L = EvlLKF(obj, theta)
            ndim = obj.Ndim;
            c = zeros(ndim, 2);
            c(1,1) = 3.5; c(1,2) = -3.5;
            r = 2;
            omega = 0.1;
            X = -((sqrt(sum((theta - c(:,1)).^2,1)) - r).^2) / (2 * omega^2);
            Y = -((sqrt(sum((theta - c(:,2)).^2,1)) - r).^2) / (2 * omega^2);
            XY = max([X; Y], [], 1);
            L = log(1 / sqrt(2 * pi * omega^2)) + XY + log(exp(X - XY) + exp(Y - XY));
        end
    end
end