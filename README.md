# SuS Bayesian — unified SeMCS (matlab)

This update retains the tested numerical algorithms and uses the common SeMCS engine and common MCMC implementation. Method-specific state and policies remain in the method classes; state evolves in G. The common core is bundled byte for byte so installation does not depend on an unreleased reliability repository.

## Qualified entrypoint

The qualified profile is `closed_ties_unit_fallback_v1`, with estimator `source_native_evidence`. The existing drivers and optional policies remain available for compatibility; use `run_qualified` for the settings in this table. No default is silently presented as paper-qualified.

```matlab
addpath(pwd); addpath(genpath(pwd));
setup_semcs();
settings = run_qualified('sus_bayes_eggbox2', true); % construct/inspect only
% result = run_qualified('sus_bayes_eggbox2'); % one actual run
```

## Source-paper comparison

Each row uses its declared independent repeat count. Elog is mean Ncal × [SD(log Z)/|printed source log Z reference|]², with sample SD (ddof=1). R=Elog/source Elog; the fixed acceptance limit is R≤1.10.

| Case | Per-level N | Repeats | Mean log Z | SD log Z | Mean Ncal | Elog | Paper Elog | R |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sus_bayes_eggbox2 | 3340 | 100 | 235.76401 | 0.29862583 | 17869.27 | 0.028645273 | 0.03724 | 0.76920712 |
| sus_bayes_loggamma2 | 5975 | 100 | -8.1953275 | 0.068951694 | 23915 | 1.6950888 | 1.686384 | 1.0051618 |
| sus_bayes_loggamma5 | 7464 | 1000 | -20.473113 | 0.16744452 | 70628.255 | 4.7259013 | 4.970112 | 0.95086415 |
| sus_bayes_shells5 | 10398 | 1000 | -5.6729835 | 0.065049242 | 51998 | 6.8439213 | 6.292 | 1.0877179 |

Source papers: [DOI 10.1016/j.ymssp.2025.112969](https://doi.org/10.1016/j.ymssp.2025.112969).

## Scope and limitations

- Qualification applies to the listed profile, estimator, targets and budgets; it does not cover every optional setting.
- No global finite-sample unbiasedness or universal posterior accuracy claim is made.
- The migration bias check uses100 same-seed native/shared pairs per case, not1000 paired runs; the1000 cohorts are independent precision studies.
- Paper efficiency points are treated as fixed; their Monte Carlo uncertainty and rounding are not included. Confidence intervals do not waive the1.10 point gate.
- Existing100-run cases are listed as100, never relabelled1000. Prior100 and1000 cohorts are not pooled.
- The closed-support/identical-seed scale fallback is an explicitly declared numerical variant, not a boundary rule claimed to have been specified by the source paper. The old author_latest default remains available and is not the qualified entrypoint.

See `VALIDATION_SUMMARY.json` for exact profile settings and source-report SHA bindings; `SEMCS_CORE_MANIFEST.json` binds every bundled core file. Existing licenses are retained, with the core notice in `LICENSE_SEMCS_CORE`. No private papers, per-run research outputs or local machine paths are bundled.
