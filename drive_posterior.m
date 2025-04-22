%% Initialization and Settings
clear;

% Simulation settings
Nsam = 1000;         % Number of samples per iteration
Nrun = 5;          % Number of repeated runs
Ndim = 2;            % Dimensionality of the parameter space

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

%% Run Simulation
Ncal = zeros(Nrun,1);            % Likelihood evaluations per run
Npos = zeros(Nrun,1);            % Number of posterior samples
Ness = zeros(Nrun,1);            % Effective sample size (posterior)

for irun = 1:Nrun
    % --- Subset Simulation ---
    out_ss = ss.RunIte;                      % Run one SuS iteration
    Ncal(irun) = sum(out_ss.Ncal);          % Total calls to likelihood

    % --- Numerical Simulation ---
    out_ps = NumSim(out_ss);                % Post-processing
    [x_pos, y_pos] = out_ps.PosSim;         % Posterior samples (U-space)
    theta_pos = IntBay.LKF.Pdis.U2X(x_pos); % Convert to parameter space

    Npos(irun) = size(x_pos,2);             % Posterior sample count

    % --- Effective Sample Size Estimation ---
    [c_z_hat, c_z] = out_ps.EvlVar;         % Estimated COVs
    Ness(irun) = (c_z / c_z_hat)^2 * Npos(irun);  % ESS under weight variance
end

%% Summary Statistics
% Effective Sample Size Ratio
Ness_rat = Ness ./ Ncal;                   % ESS per likelihood evaluation
Ness_rat_mu = mean(Ness_rat);              % Mean ESS ratio across runs