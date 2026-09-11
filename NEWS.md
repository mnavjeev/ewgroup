# ewgroup 0.2.0

- Make SURE mixing invariant to common changes of measurement units, and
  accumulate smoothing displacements directly for better numerical accuracy.
- Stream the compiled and pure-R calculations. Derivatives alone no longer
  allocate quadratic workspace; full weights are allocated only on request.
- Validate nonempty, consistently dimensioned covariance inputs using relative
  symmetry/PSD tolerances. Reject nonfinite numerical intermediates explicitly.
- Add allocation guards and interrupt checks, plus independent equation,
  finite-difference, invariance, and memory regression tests.

The public R interface and the draft's estimator and tuning default are retained.
