#' Exponentially Weighted Grouped-Heterogeneity Estimator
#'
#' Estimate grouped heterogeneity from noisy preliminary cell-specific estimates.
#' The function implements the exponentially weighted estimator and SURE
#' recombination described in the accompanying paper. Users provide the naive
#' estimates and scaled covariance estimates; the function returns the final
#' adaptive estimator.
#'
#' @param beta_hat A numeric vector of length `J` for scalar effects or a
#'   numeric `J x d` matrix for vector-valued coefficients.
#' @param Sigma_hat Scaled covariance estimates. For scalar effects, provide a
#'   numeric vector of length `J`. For vector coefficients, provide either a
#'   length-`J` list of `d x d` matrices or a `d x d x J` array.
#' @param sigma2 Positive scalar scale parameter. The covariance of
#'   `beta_hat[j, ]` is interpreted as `sigma2 * Sigma_hat[[j]]`.
#' @param gamma Positive scalar tuning parameter. If `NULL`, the default is
#'   `0.2 / max_j lambda_max(Sigma_hat[[j]])`.
#' @param return_weights Logical; if `TRUE`, include the full `J x J` matrix of
#'   exponential weights in the returned object.
#' @param return_derivative Logical; if `TRUE`, include the cell-level derivative
#'   object in the returned object and diagnostics. The default `FALSE` skips
#'   materializing full Jacobian matrices when only the final estimate and SURE
#'   diagnostics are needed.
#'
#' @return An object of class `"ewgroup_fit"` with the final estimate `theta`,
#'   the exponentially weighted estimate `tilde`, the SURE mixing weight
#'   `alpha`, the tuning parameter `gamma`, and diagnostics.
#' @aliases coef.ewgroup_fit fitted.ewgroup_fit print.ewgroup_fit
#' @importFrom stats coef fitted
#' @export
#'
#' @examples
#' beta_hat <- c(-1.2, -1.0, 0.4, 0.35, 1.1)
#' Sigma_hat <- c(1.1, 0.9, 1.0, 1.2, 0.95)
#' fit <- ewgroup(beta_hat, Sigma_hat, sigma2 = 0.05)
#' coef(fit)
ewgroup <- function(beta_hat, Sigma_hat, sigma2, gamma = NULL,
                    return_weights = FALSE, return_derivative = FALSE) {
  call <- match.call()
  beta <- as_beta_matrix(beta_hat)
  B <- beta$matrix
  J <- nrow(B)
  d <- ncol(B)

  if (!is_scalar_logical(return_weights)) {
    stop("return_weights must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is_scalar_logical(return_derivative)) {
    stop("return_derivative must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(sigma2) || length(sigma2) != 1L ||
      !is.finite(sigma2) || sigma2 <= 0) {
    stop("sigma2 must be a positive scalar.", call. = FALSE)
  }

  covariances <- normalize_covariances(Sigma_hat, J = J, d = d)
  if (is.null(gamma)) {
    gamma <- ewgroup_gamma(covariances$matrices)
  } else if (!is.numeric(gamma) || length(gamma) != 1L ||
             !is.finite(gamma) || gamma <= 0) {
    stop("gamma must be a positive scalar.", call. = FALSE)
  }

  max_lambda <- max(covariances$lambda_max)
  if (gamma * max_lambda >= 1) {
    stop(
      "gamma must satisfy gamma * max_j lambda_max(Sigma_hat[[j]]) < 1.",
      call. = FALSE
    )
  }

  core <- ewgroup_core(
    B = B,
    Sigma = covariances$matrices,
    sigma2 = sigma2,
    gamma = gamma,
    return_jacobian = return_derivative,
    return_weights = return_weights
  )

  theta <- core$alpha * core$tilde + (1 - core$alpha) * B
  derivative <- if (!return_derivative) {
    NULL
  } else if (d == 1L) {
    vapply(core$jacobian, function(x) x[1L, 1L], numeric(1L))
  } else {
    core$jacobian
  }

  out <- list(
    theta = restore_beta_shape(theta, beta),
    tilde = restore_beta_shape(core$tilde, beta),
    beta_hat = restore_beta_shape(B, beta),
    alpha = core$alpha,
    alpha_unconstrained = core$alpha_unconstrained,
    gamma = gamma,
    sigma2 = sigma2,
    sure_A = core$sure_A,
    sure_D = core$sure_D,
    derivative = derivative,
    diagnostics = list(
      max_lambda = max_lambda,
      lambda_max = covariances$lambda_max,
      sure_A = core$sure_A,
      sure_D = core$sure_D,
      jacobian = if (return_derivative) core$jacobian else NULL
    ),
    call = call
  )

  if (return_weights) {
    out$weights <- core$weights
    row_labels <- rownames(B)
    if (!is.null(row_labels)) {
      dimnames(out$weights) <- list(row_labels, row_labels)
    }
  }

  class(out) <- "ewgroup_fit"
  out
}

#' Data-Driven Default for Gamma
#'
#' Compute the default tuning parameter used by [ewgroup()].
#'
#' @param Sigma_hat Scaled covariance estimates in any format accepted by
#'   [ewgroup()].
#' @param factor Positive scalar multiplier. The default `0.2` corresponds to
#'   `gamma = 0.2 / max_j lambda_max(Sigma_hat[[j]])`.
#'
#' @return A positive scalar tuning parameter.
#' @export
#'
#' @examples
#' ewgroup_gamma(c(0.9, 1.1, 1.0))
ewgroup_gamma <- function(Sigma_hat, factor = 0.2) {
  if (!is.numeric(factor) || length(factor) != 1L ||
      !is.finite(factor) || factor <= 0) {
    stop("factor must be a positive scalar.", call. = FALSE)
  }
  covariances <- normalize_covariances(Sigma_hat)
  max_lambda <- max(covariances$lambda_max)
  if (!is.finite(max_lambda) || max_lambda <= 0) {
    stop(
      "Sigma_hat must have at least one positive eigenvalue to compute gamma.",
      call. = FALSE
    )
  }
  factor / max_lambda
}

#' @export
coef.ewgroup_fit <- function(object,
                             component = c("theta", "tilde", "naive"), ...) {
  component <- match.arg(component)
  switch(
    component,
    theta = object$theta,
    tilde = object$tilde,
    naive = object$beta_hat
  )
}

#' @export
fitted.ewgroup_fit <- function(object, ...) {
  object$theta
}

#' @export
print.ewgroup_fit <- function(x, ...) {
  theta_dim <- dim(as.matrix(x$theta))
  cat("Exponentially weighted grouped estimator\n")
  cat("  cells: ", theta_dim[1L], "\n", sep = "")
  cat("  dimension: ", theta_dim[2L], "\n", sep = "")
  cat("  gamma: ", format(signif(x$gamma, 6L)), "\n", sep = "")
  cat("  alpha: ", format(signif(x$alpha, 6L)), "\n", sep = "")
  invisible(x)
}

ewgroup_core <- function(B, Sigma, sigma2, gamma, return_jacobian = TRUE,
                         return_weights = TRUE) {
  if (exists("ewgroup_core_cpp", mode = "function")) {
    return(ewgroup_core_cpp(
      B = B,
      Sigma_hat = covariance_list_to_array(Sigma),
      sigma2 = sigma2,
      gamma = gamma,
      return_jacobian = return_jacobian,
      return_weights = return_weights
    ))
  }
  ewgroup_core_r(B = B, Sigma = Sigma, sigma2 = sigma2, gamma = gamma)
}

ewgroup_core_r <- function(B, Sigma, sigma2, gamma) {
  J <- nrow(B)
  d <- ncol(B)
  I <- diag(d)

  if (covariance_list_is_diagonal(Sigma)) {
    return(ewgroup_core_diagonal(B = B, Sigma = Sigma, sigma2 = sigma2, gamma = gamma))
  }

  Omega <- vector("list", J)
  for (k in seq_len(J)) {
    Omega[[k]] <- gamma * spd_inverse(I - gamma * Sigma[[k]])
    Omega[[k]] <- symmetrize_matrix(Omega[[k]])
  }

  weights <- matrix(NA_real_, nrow = J, ncol = J)
  tilde <- matrix(NA_real_, nrow = J, ncol = d)
  jacobian <- vector("list", J)
  trace_terms <- numeric(J)

  for (j in seq_len(J)) {
    bj <- B[j, ]
    log_weights <- numeric(J)

    for (k in seq_len(J)) {
      diff <- bj - B[k, ]
      log_weights[k] <- -0.5 * as.numeric(
        crossprod(diff, Omega[[k]] %*% diff)
      ) / sigma2
    }

    log_weights <- log_weights - max(log_weights)
    w <- exp(log_weights)
    w <- w / sum(w)
    weights[j, ] <- w

    Gamma <- weighted_matrix_sum(w, Sigma) * gamma
    H <- spd_inverse(I - Gamma)
    delta <- sweep(B, 2L, bj, "-")
    N <- as.numeric(crossprod(w, delta))
    tilde[j, ] <- bj + as.numeric(H %*% N)

    scores <- matrix(0, nrow = J, ncol = d)
    for (k in seq_len(J)) {
      diff <- bj - B[k, ]
      scores[k, ] <- -as.numeric(Omega[[k]] %*% diff) / sigma2
    }
    score_bar <- as.numeric(crossprod(w, scores))

    P <- matrix(NA_real_, nrow = d, ncol = d)
    for (r in seq_len(d)) {
      h <- numeric(d)
      h[r] <- 1
      dw <- w * (scores[, r] - score_bar[r])
      dN <- as.numeric(crossprod(dw, delta)) - (1 - w[j]) * h
      dGamma <- weighted_matrix_sum(dw, Sigma) * gamma
      P[, r] <- h +
        as.numeric(H %*% dGamma %*% H %*% N) +
        as.numeric(H %*% dN)
    }

    jacobian[[j]] <- P
    trace_terms[j] <- sum(diag(Sigma[[j]] %*% (P - I)))
  }

  sure_A <- sum((tilde - B)^2)
  sure_D <- sigma2 * sum(trace_terms)
  alpha_unconstrained <- if (sure_A <= .Machine$double.eps) {
    0
  } else {
    -sure_D / sure_A
  }
  alpha <- min(1, max(0, alpha_unconstrained))

  list(
    tilde = tilde,
    alpha = alpha,
    alpha_unconstrained = alpha_unconstrained,
    sure_A = sure_A,
    sure_D = sure_D,
    weights = weights,
    jacobian = jacobian
  )
}

ewgroup_core_diagonal <- function(B, Sigma, sigma2, gamma) {
  J <- nrow(B)
  d <- ncol(B)
  Sigma_diag <- matrix(0, nrow = J, ncol = d)
  for (j in seq_len(J)) {
    Sigma_diag[j, ] <- diag(Sigma[[j]])
  }

  Omega_diag <- gamma / (1 - gamma * Sigma_diag)
  weights <- matrix(NA_real_, nrow = J, ncol = J)
  log_weights <- matrix(0, nrow = J, ncol = J)

  for (k in seq_len(J)) {
    diff <- B
    for (coordinate in seq_len(d)) {
      diff[, coordinate] <- diff[, coordinate] - B[k, coordinate]
    }
    log_weights[, k] <- -0.5 * rowSums(
      sweep(diff^2, 2L, Omega_diag[k, ], "*")
    ) / sigma2
  }

  log_weights <- sweep(log_weights, 1L, apply(log_weights, 1L, max), "-")
  weights <- exp(log_weights)
  weights <- weights / rowSums(weights)

  Gamma <- gamma * weights %*% Sigma_diag
  H <- 1 / (1 - Gamma)
  N <- weights %*% B - B
  tilde <- B + H * N

  jacobian <- vector("list", J)
  trace_terms <- numeric(J)
  delta_coordinate <- numeric(J)
  score_coordinate <- numeric(J)
  dw <- numeric(J)

  for (j in seq_len(J)) {
    w <- weights[j, ]
    P <- matrix(0, nrow = d, ncol = d)
    delta <- sweep(B, 2L, B[j, ], "-")

    for (r in seq_len(d)) {
      delta_coordinate <- delta[, r]
      score_coordinate <- Omega_diag[, r] * delta_coordinate / sigma2
      score_bar <- sum(w * score_coordinate)
      dw <- w * (score_coordinate - score_bar)

      dN <- as.numeric(crossprod(dw, delta))
      dN[r] <- dN[r] - (1 - w[j])
      dGamma <- gamma * as.numeric(crossprod(dw, Sigma_diag))

      P[, r] <- H[j, ] * (dGamma * H[j, ] * N[j, ] + dN)
      P[r, r] <- P[r, r] + 1
    }

    jacobian[[j]] <- P
    trace_terms[j] <- sum(Sigma_diag[j, ] * (diag(P) - 1))
  }

  sure_A <- sum((tilde - B)^2)
  sure_D <- sigma2 * sum(trace_terms)
  alpha_unconstrained <- if (sure_A <= .Machine$double.eps) {
    0
  } else {
    -sure_D / sure_A
  }
  alpha <- min(1, max(0, alpha_unconstrained))

  list(
    tilde = tilde,
    alpha = alpha,
    alpha_unconstrained = alpha_unconstrained,
    sure_A = sure_A,
    sure_D = sure_D,
    weights = weights,
    jacobian = jacobian
  )
}

as_beta_matrix <- function(beta_hat) {
  if (is.numeric(beta_hat) && is.null(dim(beta_hat))) {
    if (any(!is.finite(beta_hat))) {
      stop("beta_hat must contain only finite values.", call. = FALSE)
    }
    B <- matrix(as.numeric(beta_hat), ncol = 1L)
    if (!is.null(names(beta_hat))) {
      rownames(B) <- names(beta_hat)
    }
    return(list(
      matrix = B,
      vector_input = TRUE,
      names = names(beta_hat),
      dimnames = NULL
    ))
  }

  if (is.matrix(beta_hat) && is.numeric(beta_hat)) {
    if (any(!is.finite(beta_hat))) {
      stop("beta_hat must contain only finite values.", call. = FALSE)
    }
    if (nrow(beta_hat) < 1L || ncol(beta_hat) < 1L) {
      stop("beta_hat must have at least one row and one column.", call. = FALSE)
    }
    return(list(
      matrix = beta_hat,
      vector_input = FALSE,
      names = NULL,
      dimnames = dimnames(beta_hat)
    ))
  }

  stop("beta_hat must be a numeric vector or numeric matrix.", call. = FALSE)
}

restore_beta_shape <- function(x, beta) {
  if (beta$vector_input) {
    out <- as.numeric(x[, 1L])
    names(out) <- beta$names
    return(out)
  }
  dimnames(x) <- beta$dimnames
  x
}

normalize_covariances <- function(Sigma_hat, J = NULL, d = NULL,
                                  tol = sqrt(.Machine$double.eps)) {
  if (is.numeric(Sigma_hat) && is.null(dim(Sigma_hat))) {
    if (!is.null(d) && d != 1L) {
      stop(
        "Sigma_hat must be a list or d x d x J array for vector beta_hat.",
        call. = FALSE
      )
    }
    if (!is.null(J) && length(Sigma_hat) != J) {
      stop("Sigma_hat must have the same length as beta_hat.", call. = FALSE)
    }
    if (any(!is.finite(Sigma_hat))) {
      stop("Sigma_hat must contain only finite values.", call. = FALSE)
    }
    if (any(Sigma_hat < -tol)) {
      stop("Sigma_hat must be positive semidefinite.", call. = FALSE)
    }
    values <- pmax(as.numeric(Sigma_hat), 0)
    return(list(
      matrices = lapply(values, function(x) matrix(x, nrow = 1L, ncol = 1L)),
      lambda_max = values
    ))
  }

  if (is.list(Sigma_hat)) {
    if (length(Sigma_hat) < 1L) {
      stop("Sigma_hat must contain at least one covariance matrix.", call. = FALSE)
    }
    if (!is.null(J) && length(Sigma_hat) != J) {
      stop("Sigma_hat must have one covariance matrix per row of beta_hat.", call. = FALSE)
    }
    matrices <- vector("list", length(Sigma_hat))
    lambda_max <- numeric(length(Sigma_hat))
    for (j in seq_along(Sigma_hat)) {
      checked <- check_covariance_matrix(Sigma_hat[[j]], d = d, tol = tol)
      matrices[[j]] <- checked$matrix
      lambda_max[j] <- checked$lambda_max
    }
    return(list(matrices = matrices, lambda_max = lambda_max))
  }

  dims <- dim(Sigma_hat)
  if (is.numeric(Sigma_hat) && length(dims) == 3L) {
    if (dims[1L] != dims[2L]) {
      stop("Sigma_hat arrays must have dimensions d x d x J.", call. = FALSE)
    }
    if (!is.null(d) && dims[1L] != d) {
      stop("Sigma_hat covariance dimension must match ncol(beta_hat).", call. = FALSE)
    }
    if (!is.null(J) && dims[3L] != J) {
      stop("Sigma_hat must have one covariance matrix per row of beta_hat.", call. = FALSE)
    }
    matrices <- vector("list", dims[3L])
    lambda_max <- numeric(dims[3L])
    for (j in seq_len(dims[3L])) {
      checked <- check_covariance_matrix(
        matrix(Sigma_hat[, , j], nrow = dims[1L], ncol = dims[2L]),
        d = d,
        tol = tol
      )
      matrices[[j]] <- checked$matrix
      lambda_max[j] <- checked$lambda_max
    }
    return(list(matrices = matrices, lambda_max = lambda_max))
  }

  stop(
    "Sigma_hat must be a numeric vector, a list of matrices, or a d x d x J array.",
    call. = FALSE
  )
}

check_covariance_matrix <- function(x, d = NULL, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Each covariance estimate must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(x) != ncol(x)) {
    stop("Each covariance estimate must be square.", call. = FALSE)
  }
  if (!is.null(d) && nrow(x) != d) {
    stop("Covariance dimensions must match ncol(beta_hat).", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("Sigma_hat must contain only finite values.", call. = FALSE)
  }

  scale <- max(1, abs(x))
  if (max(abs(x - t(x))) > tol * scale) {
    stop("Each covariance estimate must be symmetric.", call. = FALSE)
  }

  x <- symmetrize_matrix(x)
  eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) < -tol * scale) {
    stop("Each covariance estimate must be positive semidefinite.", call. = FALSE)
  }

  list(matrix = x, lambda_max = max(eigenvalues))
}

weighted_matrix_sum <- function(weights, matrices) {
  out <- matrices[[1L]] * weights[1L]
  if (length(matrices) > 1L) {
    for (j in 2L:length(matrices)) {
      out <- out + weights[j] * matrices[[j]]
    }
  }
  out
}

covariance_list_to_array <- function(Sigma) {
  d <- nrow(Sigma[[1L]])
  J <- length(Sigma)
  out <- array(0, dim = c(d, d, J))
  for (j in seq_len(J)) {
    out[, , j] <- Sigma[[j]]
  }
  out
}

symmetrize_matrix <- function(x) {
  0.5 * (x + t(x))
}

covariance_list_is_diagonal <- function(Sigma, tol = sqrt(.Machine$double.eps)) {
  d <- nrow(Sigma[[1L]])
  if (d == 1L) {
    return(TRUE)
  }
  off_diagonal <- row(Sigma[[1L]]) != col(Sigma[[1L]])
  for (j in seq_along(Sigma)) {
    scale <- max(1, abs(Sigma[[j]]))
    if (max(abs(Sigma[[j]][off_diagonal])) > tol * scale) {
      return(FALSE)
    }
  }
  TRUE
}

spd_inverse <- function(x) {
  chol2inv(chol(symmetrize_matrix(x)))
}

is_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}
