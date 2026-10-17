# Subset Simulation for Bayesian Updating

MATLAB implementation of subset simulation (SuS) for Bayesian evidence estimation and posterior sampling, using a Nataf transformation to work in standard normal space. This code accompanies the paper:

> **From Likelihood to Limit State: A Reliability-Inspired Framework for Bayesian Evidence Estimation and High-dimensional Sampling.**

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Benchmark likelihoods](#benchmark-likelihoods)
- [Outputs](#outputs)
- [Citation](#citation)
- [License](#license)

## Overview

The Bayesian updating problem is reformulated in the spirit of a structural reliability problem: the likelihood plays the role of a limit-state function, and subset simulation generates samples from a sequence of nested likelihood levels, from the prior toward the high-likelihood region. Samples at each level are generated with adaptive conditional sampling (aCS) MCMC. The samples from all levels are then post-processed to

- estimate the log-evidence, together with its coefficient of variation (CoV), with a correction for the correlation of Markov chain samples, and
- obtain posterior samples by resampling.

## Requirements

- MATLAB
- Statistics and Machine Learning Toolbox (`makedist` and probability distribution objects)

## Repository structure

```
.
├── drive_evidence.m                     % driver: log-evidence estimation over repeated runs
├── drive_posterior.m                    % driver: posterior sampling and effective sample size
├── SuS_Bay_Nataf.m                      % core: subset simulation for Bayesian updating (intermediate levels, seeds, convergence)
├── main.m                               % sequential sampling loop (sample generation, level update, seed selection)
├── NumSim.m                             % post-processing: evidence, posterior resampling, CoV estimates
├── Nataf.m                              % Nataf transformation between physical and standard normal space
├── MCMC-samplers/
│   └── aCS.m                            % adaptive conditional sampling
├── Likelihood-functions/                % benchmark likelihoods
├── support functions/                   % thresholds, weights, CoV correction, plotting and other helpers
├── Operations on Logarithmic values/    % log-domain arithmetic (logsum, logplus, logmean, ...)
├── LICENSE
└── README.md
```

## Quick start

1. Add the repository and all its subfolders to the MATLAB path, for example from the repository root:

   ```matlab
   >> addpath(genpath(pwd))
   ```

2. Run one of the drivers:

   ```matlab
   >> drive_evidence    % Bayesian evidence estimation
   >> drive_posterior   % posterior sampling
   ```

Both drivers use 1000 samples per level. `drive_evidence` runs 500 repetitions with `Ndim = 5`, and `drive_posterior` runs 5 repetitions with `Ndim = 2`. To change the problem, edit `Ndim` and the `Tfun = ...` line in the *Define Likelihood Function* section of the driver.

## Benchmark likelihoods

| File | Description | Dimension | Reference log-evidence |
|---|---|---|---|
| `LKF_nD_GaussianShells.m` | two Gaussian shells (multimodal), uniform prior on [−6, 6]ⁿ | n | −1.75 (n = 2), −5.67 (n = 5), −14.59 (n = 10), −36.09 (n = 20), −60.13 (n = 30) |
| `LKF_nD_GaussianLogGamma.m` | Gaussian + LogGamma mixture (skewed), uniform prior on [−30, 30]ⁿ | n | −8.19 (n = 2), −20.47 (n = 5), −40.94 (n = 10), −81.89 (n = 20), −122.82 (n = 30) |
| `LKF_2D_Eggbox.m` | eggbox (multimodal), uniform prior on [0, 10π]² | 2 | 235.856 |

The reference values are listed in the header of each file. A custom problem can be defined by writing a class with the same interface: a prior `Pdis` (a `Nataf` object), the dimension `Ndim`, and a method `EvlLKF(theta)` returning the log-likelihood of samples `theta` (`Ndim × Nsam`).

## Outputs

During each run, the progress of every level is printed (sample size, number of chains, samples per chain, likelihood evaluations). After all runs, the results are left in the workspace:

**`drive_evidence`**

| Variable | Meaning |
|---|---|
| `logz` | estimated log-evidence of each run |
| `logz_mu`, `logz_sig` | mean and standard deviation of the log-evidence |
| `cv_logz` | coefficient of variation of the log-evidence |
| `Ncal_mu` | average number of likelihood evaluations per run |

**`drive_posterior`**

| Variable | Meaning |
|---|---|
| `theta_pos` | posterior samples in the parameter space (last run) |
| `Npos` | number of posterior samples of each run |
| `Ness` | effective sample size of each run |
| `Ness_rat_mu` | average effective sample size per likelihood evaluation |

## Citation

If you use this code in your research, please cite:

> From Likelihood to Limit State: A Reliability-Inspired Framework for Bayesian Evidence Estimation and High-dimensional Sampling.

## License

This project is licensed under the MIT License. See [LICENSE](./LICENSE) for details.
