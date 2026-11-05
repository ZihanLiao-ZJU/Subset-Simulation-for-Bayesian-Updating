# Subset Simulation for Bayesian Updating

MATLAB toolbox for estimating Bayesian evidence and drawing posterior samples.

## Quick start

Open MATLAB in this repository. The common SeMCS core is included; no separate core checkout is needed.

```matlab
addpath(genpath(pwd));
setup_semcs();
rng(1);
result = run_example('sus_bayes_loggamma2');
logZ = result.logZ;
Ncal = result.Ncal;
```

To inspect the example settings before sampling:

```matlab
settings = run_example('sus_bayes_loggamma2', true);
```

Available example names: `sus_bayes_eggbox2`, `sus_bayes_loggamma2`, `sus_bayes_loggamma5`, `sus_bayes_shells5`.

## Architecture

- `+semcs/`: common sequential engine and MCMC implementations.
- `MCMC-samplers/`: public sampler classes backed by the common core.
- `Likelihood-functions/`: example likelihoods and prior distributions.
- `support functions/`: numerical integration, weights and supporting operations.
- `SuS_Bay_Nataf.m`: intermediate distributions, seed selection and stopping rules.
- `G`: evolving method parameters and sampler settings.

The example runner uses closed threshold support, aCS and native `NumSim` evidence integration. Identical seed vectors temporarily use a unit proposal scale.

The existing `drive_evidence.m` and `drive_posterior.m` scripts remain available.

## Reference

[Source paper](https://doi.org/10.1016/j.ymssp.2025.112969).

MIT license: see `LICENSE`; the bundled common core retains its notice in `LICENSE_SEMCS_CORE`.
