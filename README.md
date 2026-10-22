# Subset Simulation for Bayesian Updating

MATLAB implementation of subset simulation (SuS) for Bayesian evidence estimation and posterior sampling, using a Nataf transformation to work in standard normal space. This code accompanies the paper:

> **From likelihood to limit state: A reliability-inspired framework for Bayesian evidence estimation and high-dimensional sampling.**
> *Mechanical Systems and Signal Processing*, 237 (2025) 112969. doi: [10.1016/j.ymssp.2025.112969](https://doi.org/10.1016/j.ymssp.2025.112969)

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Repository structure](#repository-structure)
- [Quick start](#quick-start)
- [Benchmark likelihoods](#benchmark-likelihoods)
- [Custom likelihood](#custom-likelihood)
- [Outputs](#outputs)
- [Citation](#citation)
- [License](#license)

## Overview

Given observed data $D$ and uncertain parameters $\theta$, Bayesian updating follows Bayes' theorem:

$$
p(\theta \mid D)
= \frac{p(D \mid \theta)\,p(\theta)}{p(D)},
\qquad
Z = p(D) = \int p(D \mid \theta)\,p(\theta)\,\mathrm{d}\theta.
$$

Here, $p(\theta)$ is the **prior distribution**, $p(D \mid \theta)$ is the **likelihood**, $p(\theta \mid D)$ is the **posterior distribution**, and $Z = p(D)$ is the **Bayesian evidence (marginal likelihood)** that normalizes the posterior.

The toolbox provides two main outputs:

1. an **evidence (marginal likelihood) estimate** $\widehat{Z}$, computed and reported in the natural-log domain as $\log \widehat{Z}$ for numerical stability;
2. **posterior distribution samples**, which give a Monte Carlo approximation to $p(\theta \mid D)$.

The evidence is an estimate, not an exact integral. Posterior samples are obtained by resampling MCMC-generated samples; they may be repeated or dependent and should not be treated as exact independent draws.

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

The reference values are listed in the header of each file. To solve your own problem, see [Custom likelihood](#custom-likelihood).

## Custom likelihood

Any user-defined Bayesian updating problem can be solved by writing a likelihood class with the same interface as the files in `Likelihood-functions/` and passing it to the driver in place of the benchmark. The algorithm only accesses the following members:

| Member | Type | Description |
|---|---|---|
| `Ndim` | scalar | number of uncertain parameters |
| `Pdis` | `Nataf` object | prior distribution, built as `Nataf(Marg_X, Rho_X)` from a column array of MATLAB distribution objects (`makedist`) and their correlation matrix |
| `EvlLKF(obj, theta)` | method | log-likelihood of the samples `theta` (`Ndim × Nsam`, in the physical parameter space); must return a `1 × Nsam` row vector |

The sampler works in standard normal space and maps samples to the parameter space through `Pdis.U2X` before calling `EvlLKF`, so the likelihood is written directly in terms of the physical parameters. `EvlLKF` is called with many samples at once and must return the **logarithm** of the likelihood. Retain the likelihood normalization constants when estimating absolute evidence: omitting a parameter-independent factor leaves the normalized posterior unchanged but rescales the evidence.

Minimal example: two parameters with independent standard normal priors and a Gaussian likelihood for two noisy observations. Save it as `MyLikelihood.m` in a folder on the MATLAB path (e.g. `Likelihood-functions/`):

```matlab
classdef MyLikelihood
    properties
        Ndim = 2;
        Pdis = Nataf(repmat(makedist("Normal","mu",0,"sigma",1),2,1), eye(2));
    end

    methods
        function logL = EvlLKF(~, theta)
            % theta : [Ndim x Nsam] parameter samples
            % logL  : [1 x Nsam] log-likelihood
            y_obs = [1.0; -0.5];   % observed data
            sigma = 0.2;           % measurement noise standard deviation
            logL = -0.5*sum((theta - y_obs).^2, 1)/sigma^2 - numel(y_obs)*log(sqrt(2*pi)*sigma);
        end
    end
end
```

Then, in `drive_evidence.m` or `drive_posterior.m`, replace the benchmark with your class:

```matlab
Tfun = MyLikelihood;
```

The rest of the driver is unchanged.

## Outputs

The post-processing object `out_ps = NumSim(out_ss)` provides both outputs from a single subset simulation run:

- **Evidence:** `out_ps.Z` is $\log \widehat{Z}$, the natural logarithm of the evidence estimate, despite the property name `Z`. `exp(out_ps.Z)` converts it to $\widehat{Z}$ when representable; keeping the log value avoids underflow or overflow.
- **Posterior samples:** `[x_pos, y_pos] = out_ps.PosSim()` returns samples in standard normal space. `theta_pos = Tfun.Pdis.U2X(x_pos)` maps them to the physical parameter space, with one posterior sample per column.

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

> Liao, Z., Li, B., & Wan, H.-P. (2025). From likelihood to limit state: A reliability-inspired framework for Bayesian evidence estimation and high-dimensional sampling. *Mechanical Systems and Signal Processing*, 237, 112969. https://doi.org/10.1016/j.ymssp.2025.112969

```bibtex
@article{liao2025likelihood,
  title   = {From likelihood to limit state: A reliability-inspired framework for {B}ayesian evidence estimation and high-dimensional sampling},
  author  = {Liao, Zihan and Li, Binbin and Wan, Hua-Ping},
  journal = {Mechanical Systems and Signal Processing},
  volume  = {237},
  pages   = {112969},
  year    = {2025},
  doi     = {10.1016/j.ymssp.2025.112969}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](./LICENSE) for details.
