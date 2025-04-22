%% Initialization and Settings
clear;

% Simulation settings
Nsam = 1000;         % Number of samples per iteration
Nrun = 500;          % Number of repeated runs
Ndim = 5;            % Dimensionality of the parameter space

%% Define Likelihood Function (Choose One)
% ----------------------------------------------------------------------
Tfun = LKF_nD_GaussianShells(Ndim);    % Multimodal shell-shaped posterior
% Tfun = LKF_2D_Eggbox();                % Eggbox test function (2D only)
% Tfun = LKF_nD_GaussianLogGamma(Ndim);  % Skewed likelihood test case
% ----------------------------------------------------------------------

%% Construct Sampler and Subset Simulation Framework
% ----------------------------------------------------------------------
IntBay = SuS_Bay_Nataf(Tfun);    % Define Bayesian model with Nataf transformation
sampler = aCS(IntBay);           % Define adaptive MCMC sampler (e.g., aCS)
ss = main(sampler);              % Main Subset Simulation object
ss.NumSam = Nsam;                % Set sample size per level
% ----------------------------------------------------------------------

%% Run Simulations
logz = zeros(Nrun,1);             % Log-evidence per run
Ncal = zeros(Nrun,1);             % Likelihood evaluations per run

for irun = 1:Nrun
    % --- Subset Simulation ---
    out_ss = ss.RunIte;                    % Run one SuS iteration
    Ncal(irun) = sum(out_ss.Ncal);         % Total calls to likelihood
    % --- Numerical Simulation ---
    out_ps = NumSim(out_ss);               % Post-processing
    logz(irun) = out_ps.Z;                 % Store estimated log-evidence
end

%% Summary Statistics
logz_mu = mean(logz);                  % Mean log-evidence
logz_sig = std(logz);                  % Std of log-evidence
cv_logz = logz_sig / logz_mu;          % Coefficient of variation
Ncal_mu = mean(Ncal);                  % Mean likelihood calls per run
