# ewgroup

`ewgroup` implements the exponentially weighted grouped-heterogeneity estimator
from the associated paper. It takes naive cell-specific estimates as input,
uses scaled covariance estimates to form exponential weights, and combines the
weighted estimator with the naive estimator using a feasible SURE rule.

## Example

```r
library(ewgroup)

beta_hat <- c(-1.2, -1.0, 0.4, 0.35, 1.1)
Sigma_hat <- c(1.1, 0.9, 1.0, 1.2, 0.95)

fit <- ewgroup(beta_hat, Sigma_hat, sigma2 = 0.05)
coef(fit)
```

For vector-valued cell coefficients, pass a `J x d` matrix of naive estimates
and either a list of covariance matrices or a `d x d x J` array:

```r
B <- rbind(
  c(-1.0, 0.2),
  c(-0.9, 0.1),
  c(0.7, -0.4)
)
Sigma <- replicate(3, diag(c(1, 1.5)), simplify = "array")

fit <- ewgroup(B, Sigma, sigma2 = 0.1)
fitted(fit)
```

If `gamma` is omitted, the package uses
`0.2 / (d * max_j lambda_max(Sigma_hat_j))`, where `d` is the coefficient
dimension.

By default, `ewgroup()` skips materializing the full cell-level derivative
object. Set `return_derivative = TRUE` when those derivative matrices are
needed.


If you have estimates and their standard errors, pass their squared standard
errors as covariances with `sigma2 = 1`:

```r
se_hat <- c(0.2, 0.25, 0.18, 0.3, 0.22)
fit <- ewgroup(beta_hat, se_hat^2, sigma2 = 1)
```

For large jobs, keep `return_weights = FALSE` (the default). Both the compiled
core and pure R fallback process one target cell at a time, using memory
proportional to `J * d^2`, including when `return_derivative = TRUE`. Computation
still grows quadratically in `J`; requesting all weights adds about `8 * J^2`
bytes (20 GB for 50,000 cells). The compiled implementation supports R's usual
interrupt command during long calculations.

Covariance matrices must be nonempty and share one dimension. Checks for
symmetry and positive semidefiniteness use `sqrt(.Machine$double.eps)` relative
to each matrix's largest absolute entry. Tiny negative eigenvalues attributable
to roundoff in full matrices are projected to zero; negative scalar or diagonal
variances are rejected. Nonzero off-diagonal entries are always retained.

The SURE calculation uses the smoothing displacement directly and rescales its
norm before choosing `alpha`. Changing common measurement units or adding a
common shift therefore preserves the smoothing decision up to rounding. With
exactly zero displacement, `alpha = 0`. Nonfinite intermediate results produce
an error asking for rescaled inputs; extremely small `sure_A` and `sure_D`
may underflow to zero while their scale-aware ratio remains usable.
