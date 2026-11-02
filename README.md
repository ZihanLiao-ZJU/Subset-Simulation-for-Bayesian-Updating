# Subset Simulation for Bayesian Updating — MATLAB

## Overview
Estimates log Bayesian evidence and obtains weighted/resampled posterior samples using subset simulation in Nataf standard-normal space. The default sampler is aCS. Optional kernels do not mean every combination was evaluated in the cited paper.

## Structure
- `drive_evidence.m` — repeated-run example driver.
- `drive_posterior.m` — posterior samples and weight-ESS diagnostics.
- `main.m` — sequential sampling loop.
- `SuS_Bay_Nataf.m` — intermediate distributions, seeds and stopping conditions.
- `NumSim.m` — evidence and posterior post-processing.
- `Likelihood-functions/` — bundled likelihood/prior models.
- `MCMC-samplers/` — four kernels and their shared `MCMCBase` protocol.

## Quick start
Run from this repository's root. Each checkout contains its own helpers and does not require a sibling repository.

```matlab
addpath(genpath(pwd));
rng(1, 'twister');
target = LKF_nD_GaussianShells(2);
simulation = main(aCS(SuS_Bay_Nataf(target)));
simulation.NumSam = 1000;
out = simulation.RunIte();
result = NumSim(out);
logZ = result.Z;
Ncal = sum(out.Ncal);
terminal_converged = out.intBay.FlgCvg;
```

Inspect both the returned loop flag and `intBay.FlgCvg` (the intermediate model's terminal condition). A zero-seed exit can precede the terminal check; a finite evidence estimate alone is not proof of convergence. Likelihood calls are reported by `Ncal`, not by the number of posterior resamples.

## Examples and samplers
| Example | Model / dimension | Selector |
|---|---|---|
| Eggbox | Two-dimensional multimodal likelihood | `LKF_2D_Eggbox()` |
| GaussianShells | Dimension supplied by the caller | `LKF_nD_GaussianShells(Ndim)` |
| GaussianLogGamma | Dimension supplied by the caller, at least 2 | `LKF_nD_GaussianLogGamma(Ndim)` |

The existing drivers use GaussianShells and 1000 samples per level: `drive_evidence` uses dimension 5 and 500 repetitions; `drive_posterior` uses dimension 2 and 5 repetitions. The default conditional probability is 0.1. The quick start is a small usage example, not a paper-reproduction experiment.

Available samplers are `CW_aMH`, `aCS`, `HMC`, and `ES`. Set `SamplerClass = @aCS` (or another listed kernel) in the driver, or replace the constructor in the example above.

All kernels implement `SamGen(Nchains, Nsteps)` or `SamGen(U_seed, Y_seed, Nsteps)` and return samples, packed outputs, evaluation count and updated Bayesian state. Arrays use `Ndim × Nchains × Nsteps`; packed output rows are the intermediate log-likelihood factor, Gaussian log density and cached model outputs. Sampler dimensions are model-defined, not hard-coded to two.

## Dependencies
MATLAB with Statistics and Machine Learning Toolbox (`makedist`, distribution objects and resampling). Focused MATLAB-source checks also used GNU Octave 9.4 with private compatibility adapters; that does not establish native MATLAB/toolbox equivalence or make those adapters a public dependency.

## Method & validation notes
- **Implementation:** Intermediate thresholds, evidence/posterior formulas and model definitions are separate from the MCMC kernels. Numeric threshold vectors are stored in `G.parameters` (Python `G["parameters"]`); adaptation state occupies separate fields. The aCS `lambda` is a sampler scale, not the SeMIS intermediate-distribution threshold.
- **Sampler policy:** aCS uses initial coordinate scale 1 and initial lambda 0.6, passing the final lambda to the next level. ES uses coordinatewise initial angles in `[0, pi]`. HMC uses analytic Gaussian-reference rotations with likelihood acceptance in Nataf space, not generic physical-space gradient integration.
- **Outputs:** The primary evidence output is `NumSim.Z` (natural log). `EvlVar` estimates evidence uncertainty; the posterior driver computes evidence-correlation-adjusted weight ESS, not an ESS for every posterior coordinate. Resampled posterior draws may repeat or remain dependent.
- **Validation scope:** Interface, same-seed trajectory and evaluation-accounting checks are not proof of agreement with all paper statistics. Historical aCS uses seed-based scales/reset adaptation, so the current policy can change termination and results. Stopping rules have not been relaxed to hide unsuccessful runs.
- **History:** Historical MATLAB source and current sampler-state wiring are audited separately. This documentation cleanup does not change the estimator, examples or driver defaults.


- **Correlation scope:** Bundled examples use independent priors. The inherited MATLAB nonstandard-correlated-marginal branch references `GauIntPot`, which is not bundled; that branch is not independently runnable without the original helper. Python has its own quadrature implementation, but preserves the original positive-only search interval rather than extending negative-correlation support.

## Reference & license
[Paper DOI: 10.1016/j.ymssp.2025.112969](https://doi.org/10.1016/j.ymssp.2025.112969)
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

Licensed under the [MIT License](LICENSE); the original copyright notice is retained.
