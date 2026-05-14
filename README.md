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
`0.2 / max_j lambda_max(Sigma_hat_j)`.
