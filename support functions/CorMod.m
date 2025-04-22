function [c_f_hat, c_f, gamma] = CorMod(f)
% Estimate the coefficient of variation (COV) for weighted samples f,
% and apply correction for autocorrelation in Markov chains.
%
% INPUT:
%   f         : log-weight matrix [Nc × Ns], where Nc is number of chains and Ns samples per chain
%
% OUTPUT:
%   c_f       : Standard coefficient of variation of exp(f)
%   c_f_hat   : Corrected coefficient of variation with autocorrelation correction
%   gamma     : Correction factor (effective sample size factor)

% Dimensions
[Nc, Ns] = size(f);      % Nc: number of chains, Ns: samples per chain
N = Nc * Ns;             % total number of samples

% =========================================================================
% Step 1: Compute raw coefficient of variation c_f
% =========================================================================

% Mean in log-domain
f_bar = logmean(f, "all");   % log(mean(exp(f)))

% Log-deviation from mean
df = logminus(f, f_bar);     % log(exp(f) - mean(exp(f)))

% Sample variance in log-domain: var(exp(f))
f_var = logsum(2 * real(df), "all") - log(N - 1);  % log(var)

% COV computation: c_f = std(exp(f)) / mean(exp(f))
if isinf(f_bar) || isnan(f_bar) || isnan(f_var)
    c_f = 0;
else
    c_f = exp((f_var - log(N)) / 2 - f_bar);
end

% =========================================================================
% Step 2: Compute corrected coefficient of variation c_f_hat
% =========================================================================

if c_f < eps || Ns == 1
    % No correction needed (Ns = 1 or nearly zero variance)
    gamma = sqrt(Ns);
else
    % --------------------------------------------
    % Lag-1 autocorrelation estimation (within-chain)
    % --------------------------------------------
    f1 = f(:, 1:Ns-1); f2 = f(:, 2:Ns);
    f1f2_mu = logmean(f1 + f2, 1);                       % log(mean(exp(f1 + f2)))
    df1f2_mu = logmean(logminus(f1f2_mu, 2 * f_bar), "all");
    rhos = min(abs(exp(df1f2_mu - f_var)), 1);           % lag-1 autocorrelation (bounded by 1)

    % Correction factor gamma due to autocorrelation
    if rhos == 1
        gamma = Ns;
    else
        gamma = 1 + 2 * rhos / (1 - rhos)^2 * ...
            (1 - rhos - 1 / Ns * (1 - rhos^Ns));
    end
end

% Final corrected coefficient of variation
c_f_hat = sqrt(gamma) * c_f;
end
