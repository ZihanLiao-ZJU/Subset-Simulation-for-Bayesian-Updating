classdef LKF_2D_Eggbox
    % Likelihood Function: 2D Eggbox-shaped multimodal surface
    %
    % Prior: Uniform(0, 10π) in each dimension
    % Analytical log-evidence: logZ ≈ 235.856
    %
    % Reference:
    %   Feroz et al. (2009), MNRAS, 398(4):1601–1614

    properties(Constant)
        Ndim = 2;
        Pdis = Nataf(repmat(makedist("Uniform", "lower", 0, "upper", 10*pi), 2, 1), eye(2));
    end

    methods
        function L = EvlLKF(~, theta)
            L = (2 + cos(theta(1,:)/2) .* cos(theta(2,:)/2)).^5;
        end
    end
end
