# Registered settings

{
  "profile": "closed_ties_unit_fallback_v1",
  "estimator": "source_native_evidence",
  "sampler": {
    "name": "aCS",
    "Sigma0Mode": "seed_std",
    "CarryLambda": false,
    "Lambda": 0.6,
    "AcpOpt": 0.44,
    "NumIte": 5,
    "SedSav": false,
    "identical_selected_vectors_unit_fallback": true
  }
}

Only the explicit run_qualified entrypoint selects these registered settings. Constructor defaults and legacy drivers retain their previous behavior.
