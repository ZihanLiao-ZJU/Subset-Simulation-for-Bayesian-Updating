function [rho_fi, gamma] = CorMod_cross(F, I)
% Estimate correlation and correction factor between two correlated functions F and I
% in log-domain, used in cross-variance estimation in subset simulation
%
% INPUT:
%   F : log-weight vector/matrix of function f(x), size [Nc × Ns]
%   I : log-weight vector/matrix of indicator or importance function i(x), size [Nc × Ns]
%
% OUTPUT:
%   rho_fi : Estimated marginal correlation coefficient between F and I (bounded by 1)
%   gamma  : Correction factor accounting for intra-chain autocorrelation
%
% NOTE:
%   - Log-domain operations are used to ensure numerical stability
%   - `logmean`, `logsum`, `logminus` must be defined externally for log-safe computation

% -----------------------------
% Preparation
[Nc, Ns] = size(F);      % Nc: number of chains, Ns: samples per chain
N = Nc * Ns;             % total number of samples

% -----------------------------
% Marginal means (log-domain)
h = logmean(F, "all");   % log-mean of F
p = logmean(I, "all");   % log-mean of I

% -----------------------------
% Log-deviation from means
dF = logminus(F, h);     % log(F) - log(mean_F)
dI = logminus(I, p);     % log(I) - log(mean_I)

% -----------------------------
% Sample variances (log-domain)
F_var = logsum(2 * real(dF), "all") - log(N - 1);  % log(var_F)
I_var = logsum(2 * real(dI), "all") - log(N - 1);  % log(var_I)

% -----------------------------
% Marginal correlation estimate (bounded by 1)
FI_mu = logmean(F + I, 1);                             % log(mean(exp(F + I)))
dFI_mu = logmean(logminus(FI_mu, h + p), "all");       % log(mean((F-h)(I-p)))
rho_fi = min(abs(exp(dFI_mu - (F_var + I_var) / 2)), 1); % bounded by 1

% -----------------------------
% Serial correlation estimate within Markov chains (optional for Ns > 1)
if Ns == 1
    rhos_fi = 0;
else
    F1 = F(:, 1:Ns-1); F2 = F(:, 2:Ns);
    I1 = I(:, 1:Ns-1); I2 = I(:, 2:Ns);

    F1I2_mu = logmean(F1 + I2, 1);
    I1F2_mu = logmean(I1 + F2, 1);
    dF1I2_mu = logmean(logminus(F1I2_mu, h + p), "all");
    dI1F2_mu = logmean(logminus(I1F2_mu, h + p), "all");

    rhos1 = min(abs(exp(dF1I2_mu - (F_var + I_var) / 2)), 1);
    rhos2 = min(abs(exp(dI1F2_mu - (F_var + I_var) / 2)), 1);
    rhos_fi = (rhos1 + rhos2) / 2;
end

% -----------------------------
% Correction factor gamma for effective sample size computation
if rhos_fi == 1
    gamma = Ns;
else
    gamma = 1 + 2 * rhos_fi / (1 - rhos_fi)^2 * (1 - rhos_fi - (1/Ns) * (1 - rhos_fi^Ns));
end
end
