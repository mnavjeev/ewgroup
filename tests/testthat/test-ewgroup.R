reference_scalar <- function(y, sigma2, Sigma_hat, gamma = 0.20) {
  Sigma_hat <- pmax(Sigma_hat, .Machine$double.eps)
  if (any(gamma * Sigma_hat >= 1)) {
    stop("gamma * Sigma_hat must be below one for all cells.")
  }

  Omega <- gamma / (1 - gamma * Sigma_hat)
  diff <- outer(y, y, "-")
  log_weights <- -0.5 * sweep(diff^2, 2, Omega, "*") / sigma2
  log_weights <- sweep(log_weights, 1, apply(log_weights, 1, max), "-")
  weights <- exp(log_weights)
  weights <- weights / rowSums(weights)

  Gamma <- gamma * as.numeric(weights %*% Sigma_hat)
  H <- 1 / (1 - Gamma)

  delta <- -diff
  N <- rowSums(weights * delta)
  tilde <- y + H * N

  score <- -sweep(diff, 2, Omega, "*") / sigma2
  score_bar <- rowSums(weights * score)
  d_weights <- weights * sweep(score, 1, score_bar, "-")
  d_N <- rowSums(d_weights * delta) - (1 - diag(weights))
  d_Gamma <- gamma * rowSums(sweep(d_weights, 2, Sigma_hat, "*"))
  derivative <- 1 + H^2 * d_Gamma * N + H * d_N

  sure_A <- sum((y - tilde)^2)
  sure_D <- sigma2 * sum(Sigma_hat * (derivative - 1))
  alpha_unconstrained <- if (sure_A <= .Machine$double.eps) 0 else -sure_D / sure_A
  alpha <- min(1, max(0, alpha_unconstrained))
  theta <- alpha * tilde + (1 - alpha) * y

  list(
    theta = theta,
    tilde = tilde,
    alpha = alpha,
    alpha_unconstrained = alpha_unconstrained,
    derivative = derivative,
    sure_A = sure_A,
    sure_D = sure_D,
    weights = weights
  )
}

test_that("scalar implementation matches simulation reference", {
  y <- c(a = -1.1, b = -0.95, c = 0.2, d = 0.27, e = 1.3)
  Sigma_hat <- c(0.8, 1.1, 0.9, 1.2, 1.0)
  ref <- reference_scalar(y, sigma2 = 0.05, Sigma_hat = Sigma_hat, gamma = 0.2)
  fit <- ewgroup(
    y,
    Sigma_hat,
    sigma2 = 0.05,
    gamma = 0.2,
    return_weights = TRUE,
    return_derivative = TRUE
  )

  expect_s3_class(fit, "ewgroup_fit")
  expect_equal(unname(fit$theta), unname(ref$theta), tolerance = 1e-12)
  expect_equal(unname(fit$tilde), unname(ref$tilde), tolerance = 1e-12)
  expect_equal(fit$alpha, ref$alpha, tolerance = 1e-12)
  expect_equal(fit$alpha_unconstrained, ref$alpha_unconstrained, tolerance = 1e-12)
  expect_equal(unname(fit$derivative), unname(ref$derivative), tolerance = 1e-12)
  expect_equal(fit$sure_A, ref$sure_A, tolerance = 1e-12)
  expect_equal(fit$sure_D, ref$sure_D, tolerance = 1e-12)
  expect_equal(unname(fit$weights), unname(ref$weights), tolerance = 1e-12)
  expect_equal(names(coef(fit)), names(y))
})

test_that("default gamma uses maximum covariance eigenvalue", {
  Sigma <- list(
    diag(c(1, 2)),
    matrix(c(2, 0.1, 0.1, 3), nrow = 2)
  )
  B <- rbind(c(0, 0), c(1, -1))
  expected <- 0.2 / max(vapply(
    Sigma,
    function(x) max(eigen(x, symmetric = TRUE, only.values = TRUE)$values),
    numeric(1)
  ))

  expect_equal(ewgroup_gamma(Sigma), expected)
  expect_equal(ewgroup(B, Sigma, sigma2 = 0.1)$gamma, expected)
})

test_that("list and array covariance inputs agree for vector estimates", {
  B <- rbind(
    cell1 = c(-1.0, 0.2),
    cell2 = c(-0.9, 0.1),
    cell3 = c(0.7, -0.4)
  )
  colnames(B) <- c("x1", "x2")
  Sigma <- list(
    matrix(c(1.0, 0.1, 0.1, 1.5), 2),
    matrix(c(1.2, 0.0, 0.0, 1.1), 2),
    matrix(c(0.9, 0.2, 0.2, 1.4), 2)
  )
  Sigma_array <- array(unlist(Sigma), dim = c(2, 2, 3))

  fit_list <- ewgroup(B, Sigma, sigma2 = 0.1, gamma = 0.05, return_weights = TRUE)
  fit_array <- ewgroup(B, Sigma_array, sigma2 = 0.1, gamma = 0.05)

  expect_equal(fit_array$theta, fit_list$theta, tolerance = 1e-12)
  expect_equal(fit_array$tilde, fit_list$tilde, tolerance = 1e-12)
  expect_equal(unname(rowSums(fit_list$weights)), rep(1, nrow(B)), tolerance = 1e-12)
  expect_true(fit_list$alpha >= 0)
  expect_true(fit_list$alpha <= 1)
  expect_equal(
    fit_list$theta,
    fit_list$alpha * fit_list$tilde + (1 - fit_list$alpha) * B,
    tolerance = 1e-12
  )
})

test_that("derivative materialization is skipped by default", {
  B <- rbind(
    cell1 = c(-1.0, 0.2),
    cell2 = c(-0.9, 0.1),
    cell3 = c(0.7, -0.4)
  )
  Sigma <- list(
    matrix(c(1.0, 0.1, 0.1, 1.5), 2),
    matrix(c(1.2, 0.0, 0.0, 1.1), 2),
    matrix(c(0.9, 0.2, 0.2, 1.4), 2)
  )

  fit_full <- ewgroup(
    B,
    Sigma,
    sigma2 = 0.1,
    gamma = 0.05,
    return_derivative = TRUE
  )
  fit_trace <- ewgroup(
    B,
    Sigma,
    sigma2 = 0.1,
    gamma = 0.05
  )

  expect_equal(fit_trace$theta, fit_full$theta, tolerance = 1e-12)
  expect_equal(fit_trace$tilde, fit_full$tilde, tolerance = 1e-12)
  expect_equal(fit_trace$alpha, fit_full$alpha, tolerance = 1e-12)
  expect_equal(fit_trace$sure_D, fit_full$sure_D, tolerance = 1e-12)
  expect_null(fit_trace$derivative)
  expect_null(fit_trace$diagnostics$jacobian)
})

test_that("one-column matrix accepts one-by-one covariance arrays", {
  y <- c(-1.1, -0.95, 0.2, 0.27, 1.3)
  Sigma_hat <- c(0.8, 1.1, 0.9, 1.2, 1.0)
  Sigma_array <- array(Sigma_hat, dim = c(1, 1, length(Sigma_hat)))

  fit_vector <- ewgroup(y, Sigma_hat, sigma2 = 0.05, gamma = 0.2)
  fit_matrix <- ewgroup(matrix(y, ncol = 1), Sigma_array, sigma2 = 0.05, gamma = 0.2)

  expect_equal(as.numeric(fit_matrix$theta[, 1]), unname(fit_vector$theta), tolerance = 1e-12)
  expect_equal(as.numeric(fit_matrix$tilde[, 1]), unname(fit_vector$tilde), tolerance = 1e-12)
  expect_equal(fit_matrix$alpha, fit_vector$alpha, tolerance = 1e-12)
  expect_equal(fit_matrix$sure_D, fit_vector$sure_D, tolerance = 1e-12)
})

test_that("analytic derivative matches finite differences", {
  B <- rbind(c(-1.0, 0.2), c(-0.8, 0.1), c(0.6, -0.3))
  Sigma <- list(
    matrix(c(1.0, 0.1, 0.1, 1.3), 2),
    matrix(c(1.1, 0.0, 0.0, 0.9), 2),
    matrix(c(0.8, 0.2, 0.2, 1.4), 2)
  )
  fit <- ewgroup(B, Sigma, sigma2 = 0.1, gamma = 0.05, return_derivative = TRUE)
  eps <- 1e-6
  j <- 2L
  h <- c(0.4, -0.7)
  B_eps <- B
  B_eps[j, ] <- B_eps[j, ] + eps * h
  fit_eps <- ewgroup(B_eps, Sigma, sigma2 = 0.1, gamma = 0.05)

  finite_diff <- (fit_eps$tilde[j, ] - fit$tilde[j, ]) / eps
  analytic <- as.numeric(fit$derivative[[j]] %*% h)
  expect_equal(finite_diff, analytic, tolerance = 1e-5)
})

test_that("multivariate analytic derivative matches coordinate finite differences", {
  B <- rbind(
    c(-1.0, 0.2, 0.4),
    c(-0.8, 0.1, 0.5),
    c(0.6, -0.3, -0.2),
    c(0.7, -0.4, -0.1)
  )
  Sigma <- list(
    matrix(c(1.0, 0.1, 0.0, 0.1, 1.3, 0.2, 0.0, 0.2, 1.1), 3),
    matrix(c(1.1, 0.0, 0.1, 0.0, 0.9, 0.0, 0.1, 0.0, 1.2), 3),
    matrix(c(0.8, 0.2, 0.0, 0.2, 1.4, 0.1, 0.0, 0.1, 1.0), 3),
    matrix(c(1.2, 0.1, 0.1, 0.1, 1.0, 0.0, 0.1, 0.0, 1.3), 3)
  )
  fit <- ewgroup(B, Sigma, sigma2 = 0.1, gamma = 0.05, return_derivative = TRUE)
  eps <- 1e-6

  for (j in seq_len(nrow(B))) {
    for (r in seq_len(ncol(B))) {
      B_eps <- B
      B_eps[j, r] <- B_eps[j, r] + eps
      fit_eps <- ewgroup(B_eps, Sigma, sigma2 = 0.1, gamma = 0.05)

      finite_diff <- (fit_eps$tilde[j, ] - fit$tilde[j, ]) / eps
      expect_equal(finite_diff, fit$derivative[[j]][, r], tolerance = 1e-5)
    }
  }
})

test_that("input validation catches incompatible inputs", {
  expect_error(ewgroup(c(1, 2), c(1), sigma2 = 0.1), "same length")
  expect_error(
    ewgroup(matrix(1:4, nrow = 2), c(1, 1), sigma2 = 0.1),
    "list or d x d x J array"
  )
  expect_error(ewgroup(c(1, 2), c(1, -0.1), sigma2 = 0.1), "positive semidefinite")
  expect_error(ewgroup(c(1, 2), c(1, 1), sigma2 = -0.1), "sigma2")
  expect_error(ewgroup(c(1, 2), c(1, 1), sigma2 = 0.1, gamma = 1), "gamma")
})
