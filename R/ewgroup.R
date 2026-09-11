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
#'   `0.2 / (d * max_j lambda_max(Sigma_hat[[j]]))`, where `d` is the
#'   coefficient dimension.
#' @param return_weights Logical; if `TRUE`, include the full `J x J` matrix of
#'   exponential weights in the returned object.
#' @param return_derivative Logical; if `TRUE`, include the cell-level derivative
#'   object in the returned object and diagnostics. The default `FALSE` skips
#'   materializing full Jacobian matrices when only the final estimate and SURE
#'   diagnostics are needed.
#'
#' @details The implementation uses memory linear in the number of cells when
#'   `return_weights = FALSE`, including when derivatives are requested. Runtime
#'   remains quadratic in the number of cells. A full weights matrix requires
#'   approximately `8 * J^2` bytes. Both the compiled and pure R implementations
#'   use this streaming calculation.
#'
#'   Covariance matrices must have a common nonzero dimension. Symmetry and
#'   positive semidefiniteness are checked relative to each matrix's largest
#'   absolute entry with tolerance `sqrt(.Machine$double.eps)`. Tolerated negative
#'   eigenvalues in full matrices are projected to zero; diagonal variances must
#'   be nonnegative. SURE uses a scaled displacement norm, with `alpha = 0` for
#'   exactly zero displacement. Nonfinite numerical intermediates produce a
#'   rescaling error. Diagnostics smaller than machine precision may underflow
#'   to zero without changing the scale-aware mixing ratio.
#'
#' @return An object of class `"ewgroup_fit"` with the final estimate `theta`,
#'   the exponentially weighted estimate `tilde`, the SURE mixing weight
#'   `alpha`, the tuning parameter `gamma`, and diagnostics.
#' @aliases coef.ewgroup_fit fitted.ewgroup_fit print.ewgroup_fit
#' The package imports `coef()` and `fitted()` from R's standard `stats`
#' package so this object works with the usual R commands for extracting
#' estimates. Those commands are standard "generic" functions: R chooses the
#' method below because the returned object has class `"ewgroup_fit"`.
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
  # Save the user's original command. The print method and debugging output can
  # then show exactly how the fit was requested.
  call <- match.call()

  # Internally the estimator always works with a numeric table:
  # one row per cell and one column per coefficient. A plain vector is treated
  # as a one-column table, and `beta` keeps enough information to restore the
  # user's original vector or matrix shape before returning results.
  beta <- as_beta_matrix(beta_hat)
  B <- beta$matrix
  J <- nrow(B)
  d <- ncol(B)

  # These two options must be single TRUE/FALSE values. This prevents mistakes
  # such as passing NA or a longer vector, where R might otherwise silently use
  # only part of the input.
  if (!is_scalar_logical(return_weights)) {
    stop("return_weights must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is_scalar_logical(return_derivative)) {
    stop("return_derivative must be TRUE or FALSE.", call. = FALSE)
  }

  # `sigma2` is the overall noise scale. It must be one positive number because
  # the cell-specific covariance information is supplied separately in
  # `Sigma_hat`.
  if (!is.numeric(sigma2) || length(sigma2) != 1L ||
      !is.finite(sigma2) || sigma2 <= 0) {
    stop("sigma2 must be a positive scalar.", call. = FALSE)
  }

  # Accept the covariance input in several user-friendly forms, then convert it
  # to one common form: a list with one covariance matrix per cell. The helper
  # also records the largest eigenvalue for each covariance matrix. An
  # eigenvalue is a number that describes how large the covariance matrix can
  # make a direction; here we only need the largest one to keep `gamma` safe.
  covariances <- normalize_covariances(Sigma_hat, J = J, d = d)

  # If the user does not choose `gamma`, use the documented default. Otherwise,
  # check that the supplied value is one positive number.
  if (is.null(gamma)) {
    gamma <- default_gamma(d, max(covariances$lambda_max), 0.2)
  } else if (!is.numeric(gamma) || length(gamma) != 1L ||
             !is.finite(gamma) || gamma <= 0) {
    stop("gamma must be a positive scalar.", call. = FALSE)
  }

  # The formulas below contain inverse matrices based on `gamma`. This check
  # rules out values of `gamma` that would make those inverses unstable or
  # impossible to compute.
  max_lambda <- max(covariances$lambda_max)
  if (gamma * max_lambda >= 1) {
    stop(
      "gamma must satisfy gamma * max_j lambda_max(Sigma_hat[[j]]) < 1.",
      call. = FALSE
    )
  }

  # Do the numerical work in the core routine. The core returns the
  # exponentially weighted estimate (`tilde`), the data-driven mixing weight
  # (`alpha`), and optional details such as the full weights and derivative
  # matrices.
  core <- ewgroup_core(
    B = B,
    Sigma = covariances$matrices,
    sigma2 = sigma2,
    gamma = gamma,
    return_jacobian = return_derivative,
    return_weights = return_weights
  )

  # The final estimate is a blend of the smoothed estimate and the original
  # cell estimate. When `alpha` is 1 the result is fully smoothed; when it is 0
  # the result is the original input.
  theta <- B + core$alpha * core$adjustment
  check_numerical_finite(theta)

  # The derivative object can be large for vector-valued estimates. Keep it out
  # of the returned object unless the user asked for it. For scalar estimates,
  # each derivative is only one number, so return a simple vector.
  derivative <- if (!return_derivative) {
    NULL
  } else if (d == 1L) {
    vapply(core$jacobian, function(x) x[1L, 1L], numeric(1L))
  } else {
    core$jacobian
  }

  # Build the result list. Values that have the same shape as `beta_hat` are
  # converted back to the user's original vector or matrix form.
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

  # The weights are a J by J table. Row j contains the weights used when
  # smoothing cell j. If the input estimates had names, reuse them as row and
  # column labels so the table is easier to read.
  if (return_weights) {
    out$weights <- core$weights
    row_labels <- rownames(B)
    if (!is.null(row_labels)) {
      dimnames(out$weights) <- list(row_labels, row_labels)
    }
  }

  # Giving the list a class lets R dispatch `print()`, `coef()`, and `fitted()`
  # to the methods defined below.
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
#'   `gamma = 0.2 / (d * max_j lambda_max(Sigma_hat[[j]]))`, where `d` is the
#'   covariance-matrix dimension.
#'
#' @return A positive scalar tuning parameter.
#' @export
#'
#' @examples
#' ewgroup_gamma(c(0.9, 1.1, 1.0))
ewgroup_gamma <- function(Sigma_hat, factor = 0.2) {
  # `factor` controls how aggressive the default smoothing can be. It is kept
  # positive so the resulting `gamma` is positive too.
  if (!is.numeric(factor) || length(factor) != 1L ||
      !is.finite(factor) || factor <= 0) {
    stop("factor must be a positive scalar.", call. = FALSE)
  }

  # Reuse the same covariance checker used by `ewgroup()`. This lets users pass
  # a vector, a list of matrices, or a 3D array here as well.
  covariances <- normalize_covariances(Sigma_hat)
  d <- nrow(covariances$matrices[[1L]])
  max_lambda <- max(covariances$lambda_max)

  default_gamma(d, max_lambda, factor)
}

default_gamma <- function(d, max_lambda, factor) {
  # If all covariance information is zero, the default formula would divide by
  # zero. In that case the data do not provide a scale for `gamma`.
  if (!is.finite(max_lambda) || max_lambda <= 0) {
    stop(
      "Sigma_hat must have at least one positive eigenvalue to compute gamma.",
      call. = FALSE
    )
  }

  # Scale by both the coefficient dimension and the largest covariance
  # eigenvalue, matching the paper's gamma proportional to 1 / d regime while
  # retaining invariance to a common rescaling of the covariance matrices.
  gamma <- (factor / d) / max_lambda
  if (!is.finite(gamma) || gamma <= 0) {
    stop("Default gamma is not representable; rescale Sigma_hat.", call. = FALSE)
  }
  gamma
}

#' @export
coef.ewgroup_fit <- function(object,
                             component = c("theta", "tilde", "naive"), ...) {
  # `match.arg()` checks the user's choice and also fills in the default choice,
  # which is "theta". It prevents misspelled component names from being silently
  # accepted.
  component <- match.arg(component)

  # Return the requested estimate:
  # theta = final blended estimate,
  # tilde = exponentially weighted estimate before final blending,
  # naive = original cell estimate supplied by the user.
  switch(
    component,
    theta = object$theta,
    tilde = object$tilde,
    naive = object$beta_hat
  )
}

#' @export
fitted.ewgroup_fit <- function(object, ...) {
  # The "fitted" value for this estimator is the final estimate.
  object$theta
}

#' @export
print.ewgroup_fit <- function(x, ...) {
  # Keep printing short: show the size of the problem and the two main tuning
  # outputs, then return the object invisibly so it can still be reused.
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
  # The C++ routine is much faster for larger inputs. `Rcpp` connects R to C++;
  # the generated wrapper is named `ewgroup_core_cpp`. If the package was loaded
  # in a way that did not compile C++, fall back to the pure R implementation
  # below so the estimator is still usable.
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

  # Pure R fallback. It is easier to read and useful during development, but it
  # is slower than the compiled C++ routine.
  ewgroup_core_r(B = B, Sigma = Sigma, sigma2 = sigma2, gamma = gamma,
                 return_jacobian = return_jacobian, return_weights = return_weights)
}

ewgroup_core_r <- function(B, Sigma, sigma2, gamma, return_jacobian = TRUE,
                           return_weights = TRUE) {
  ewgroup_core_stream(B, Sigma, sigma2, gamma, return_jacobian, return_weights,
                     diagonal = covariance_list_is_diagonal(Sigma))
}

ewgroup_core_diagonal <- function(B, Sigma, sigma2, gamma,
                                  return_jacobian = TRUE, return_weights = TRUE) {
  ewgroup_core_stream(B, Sigma, sigma2, gamma, return_jacobian, return_weights,
                     diagonal = TRUE)
}

ewgroup_core_stream <- function(B, Sigma, sigma2, gamma, return_jacobian,
                               return_weights, diagonal) {
  J <- nrow(B)
  d <- ncol(B)
  I <- diag(d)
  scaled_sigma <- lapply(Sigma, function(x) gamma * x)
  if (diagonal) {
    Sdiag <- matrix(unlist(lapply(Sigma, diag)), nrow = J, byrow = TRUE)
    Cdiag <- gamma * Sdiag
    Omega_diag <- gamma / (1 - Cdiag)
    check_numerical_finite(Omega_diag)
  } else {
    Omega <- lapply(scaled_sigma, function(x) gamma * spd_inverse(I - x))
    for (x in Omega) check_numerical_finite(x)
  }
  # Keep one row of pairwise weights and scores, regardless of diagnostics.
  # Only explicitly requested weights require a J by J matrix.
  if (return_weights && as.double(J) * J > 2^52) {
    stop("Requested weights matrix is too large.", call. = FALSE)
  }
  weights <- if (return_weights) matrix(0, J, J) else NULL
  adjustment <- tilde <- matrix(0, J, d)
  jacobian <- vector("list", if (return_jacobian) J else 0L)
  trace_terms <- derivative_trace <- numeric(J)
  for (j in seq_len(J)) {
    delta <- sweep(B, 2L, B[j, ], "-")
    check_numerical_finite(delta)
    if (diagonal) {
      scores <- Omega_diag * delta / sigma2
    } else {
      scores <- matrix(0, J, d)
      for (k in seq_len(J)) scores[k, ] <- Omega[[k]] %*% delta[k, ] / sigma2
    }
    check_numerical_finite(scores)
    log_weights <- -0.5 * rowSums(delta * scores)
    check_numerical_finite(log_weights)
    w <- exp(log_weights - max(log_weights))
    w <- w / sum(w)
    if (return_weights) weights[j, ] <- w
    N <- as.numeric(crossprod(w, delta))
    if (diagonal) {
      Hdiag <- 1 / (1 - as.numeric(crossprod(w, Cdiag)))
      s <- Hdiag * N
      sigma_s <- sweep(Cdiag, 2L, s, "*")
    } else {
      H <- spd_inverse(I - weighted_matrix_sum(w, scaled_sigma))
      s <- as.numeric(H %*% N)
      sigma_s <- matrix(0, J, d)
      for (k in seq_len(J)) sigma_s[k, ] <- scaled_sigma[[k]] %*% s
    }
    check_numerical_finite(s)
    adjustment[j, ] <- s
    tilde[j, ] <- B[j, ] + s
    check_numerical_finite(tilde[j, ])
    score_bar <- as.numeric(crossprod(w, scores))
    dw <- w * sweep(scores, 2L, score_bar, "-")
    U <- crossprod(delta + sigma_s, dw) - (1 - w[j]) * I
    U <- if (diagonal) Hdiag * U else H %*% U
    check_numerical_finite(U)
    trace_terms[j] <- sum(Sigma[[j]] * t(U))
    derivative_trace[j] <- sum(diag(Sigma[[j]])) + trace_terms[j]
    if (return_jacobian) jacobian[[j]] <- I + U
  }
  mix <- sure_mixing(adjustment, trace_terms, sigma2)
  c(list(tilde = tilde, adjustment = adjustment), mix,
    list(weights = weights, jacobian = jacobian, derivative_trace = derivative_trace))
}

check_numerical_finite <- function(x) {
  if (any(!is.finite(x))) {
    stop("Nonfinite numerical result; rescale beta_hat and its covariance.", call. = FALSE)
  }
  invisible(NULL)
}

sure_mixing <- function(adjustment, trace_terms, sigma2) {
  scale <- max(abs(adjustment))
  trace_sum <- sum(trace_terms)
  check_numerical_finite(c(scale, trace_sum))
  norm_sum <- if (scale == 0) 0 else sum((adjustment / scale)^2)
  sure_A <- if (scale == 0) 0 else exp(2 * log(scale) + log(norm_sum))
  sure_D <- sigma2 * trace_sum
  check_numerical_finite(c(sure_A, sure_D))
  # Logarithms preserve the mixing ratio when squaring very small movements
  # underflows. Exactly zero displacement is the only zero-norm special case.
  alpha_unconstrained <- if (scale == 0 || trace_sum == 0) 0 else {
    -sign(trace_sum) * exp(log(sigma2) + log(abs(trace_sum)) -
                           2 * log(scale) - log(norm_sum))
  }
  list(alpha = min(1, max(0, alpha_unconstrained)),
       alpha_unconstrained = alpha_unconstrained, sure_A = sure_A, sure_D = sure_D)
}

as_beta_matrix <- function(beta_hat) {
  # Accept a plain numeric vector for scalar estimates. R vectors do not have
  # rows and columns, so convert them to a one-column matrix for the core code.
  if (is.numeric(beta_hat) && is.null(dim(beta_hat))) {
    if (length(beta_hat) < 1L) {
      stop("beta_hat must contain at least one estimate.", call. = FALSE)
    }
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

  # Accept a numeric matrix for vector-valued estimates. Each row is one cell;
  # each column is one coefficient for that cell.
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

  # Anything else, such as a data frame or character vector, is rejected because
  # the estimator needs numeric cell estimates.
  stop("beta_hat must be a numeric vector or numeric matrix.", call. = FALSE)
}

restore_beta_shape <- function(x, beta) {
  # If the user supplied a vector, return a vector with the original names.
  if (beta$vector_input) {
    out <- as.numeric(x[, 1L])
    names(out) <- beta$names
    return(out)
  }

  # If the user supplied a matrix, keep the original row and column labels.
  dimnames(x) <- beta$dimnames
  x
}

normalize_covariances <- function(Sigma_hat, J = NULL, d = NULL,
                                  tol = sqrt(.Machine$double.eps)) {
  # Scalar estimates can use a simple vector of covariance estimates. Each
  # number becomes a one-by-one covariance matrix.
  if (is.numeric(Sigma_hat) && is.null(dim(Sigma_hat))) {
    if (length(Sigma_hat) < 1L) {
      stop("Sigma_hat must contain at least one covariance estimate.", call. = FALSE)
    }
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
    if (any(Sigma_hat < 0)) {
      stop("Sigma_hat must be positive semidefinite.", call. = FALSE)
    }
    values <- pmax(as.numeric(Sigma_hat), 0)
    return(list(
      matrices = lapply(values, function(x) matrix(x, nrow = 1L, ncol = 1L)),
      lambda_max = values
    ))
  }

  # Vector estimates may pass a list, where each list entry is one covariance
  # matrix for one cell. Check every matrix and collect the largest eigenvalue
  # from each one.
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
      if (is.null(d)) d <- nrow(checked$matrix)
      matrices[[j]] <- checked$matrix
      lambda_max[j] <- checked$lambda_max
    }
    return(list(matrices = matrices, lambda_max = lambda_max))
  }

  # A 3D array is another compact way to provide the same information: the first
  # two dimensions are the covariance matrix, and the third dimension indexes
  # the cell.
  dims <- dim(Sigma_hat)
  if (is.numeric(Sigma_hat) && length(dims) == 3L) {
    if (any(dims < 1L)) {
      stop("Sigma_hat arrays must have nonempty covariance dimensions and cells.", call. = FALSE)
    }
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
      if (is.null(d)) d <- nrow(checked$matrix)
      matrices[[j]] <- checked$matrix
      lambda_max[j] <- checked$lambda_max
    }
    return(list(matrices = matrices, lambda_max = lambda_max))
  }

  # At this point none of the accepted input forms matched.
  stop(
    "Sigma_hat must be a numeric vector, a list of matrices, or a d x d x J array.",
    call. = FALSE
  )
}

check_covariance_matrix <- function(x, d = NULL, tol = sqrt(.Machine$double.eps)) {
  # A covariance matrix must be a numeric square table. Square means it has the
  # same number of rows and columns.
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("Each covariance estimate must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop("Each covariance estimate must be nonempty.", call. = FALSE)
  }
  if (nrow(x) != ncol(x)) {
    stop("Each covariance estimate must be square.", call. = FALSE)
  }
  if (!is.null(d) && nrow(x) != d) {
    stop("Covariance dimensions must agree and match ncol(beta_hat) when supplied.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("Sigma_hat must contain only finite values.", call. = FALSE)
  }

  # Covariance matrices should be symmetric, meaning entry (row, column) equals
  # entry (column, row). The tolerance allows tiny roundoff differences.
  scale <- max(abs(x))
  if (scale == 0) return(list(matrix = x, lambda_max = 0))
  normalized <- x / scale
  if (max(abs(normalized - t(normalized))) > tol) {
    stop("Each covariance estimate must be symmetric.", call. = FALSE)
  }

  # Average the matrix with its transpose to remove harmless roundoff asymmetry,
  # then use `eigen()` to check that the matrix is positive semidefinite. Here
  # `eigen()` breaks the matrix into direction sizes; negative direction sizes
  # would mean an invalid covariance matrix.
  normalized <- symmetrize_matrix(normalized)
  if (all(normalized[row(x) != col(x)] == 0) && any(diag(x) < 0)) {
    stop("Each covariance estimate must be positive semidefinite.", call. = FALSE)
  }
  eig <- eigen(normalized, symmetric = TRUE)
  if (min(eig$values) < -tol) {
    stop("Each covariance estimate must be positive semidefinite.", call. = FALSE)
  }
  # Roundoff-size negative eigenvalues are projected to zero consistently with
  # the Stata implementation. Never silently discard a genuine off-diagonal.
  if (any(eig$values < 0)) {
    normalized <- tcrossprod(sweep(eig$vectors, 2L, sqrt(pmax(eig$values, 0)), "*"))
  }
  x <- normalized * scale
  lambda_max <- max(eig$values) * scale
  if (!is.finite(lambda_max)) {
    stop("Covariance eigenvalue is not representable; rescale Sigma_hat.", call. = FALSE)
  }
  list(matrix = x, lambda_max = lambda_max)
}

weighted_matrix_sum <- function(weights, matrices) {
  # Add matrices together after multiplying each one by its weight. This is the
  # matrix version of a weighted average.
  out <- matrices[[1L]] * weights[1L]
  if (length(matrices) > 1L) {
    for (j in 2L:length(matrices)) {
      out <- out + weights[j] * matrices[[j]]
    }
  }
  out
}

covariance_list_to_array <- function(Sigma) {
  # The C++ routine receives covariance matrices as one 3D array. This helper
  # stacks the R list of matrices into that shape.
  d <- nrow(Sigma[[1L]])
  J <- length(Sigma)
  out <- array(0, dim = c(d, d, J))
  for (j in seq_len(J)) {
    out[, , j] <- Sigma[[j]]
  }
  out
}

symmetrize_matrix <- function(x) {
  # Replace a matrix by the average of itself and its transpose. This keeps the
  # intended symmetric covariance matrix even if floating point arithmetic made
  # two mirrored entries differ by a tiny amount.
  0.5 * x + 0.5 * t(x)
}

covariance_list_is_diagonal <- function(Sigma) {
  # Check whether all off-diagonal entries are exactly zero. If so, the
  # estimator can use the faster diagonal calculation above.
  d <- nrow(Sigma[[1L]])
  if (d == 1L) {
    return(TRUE)
  }
  off_diagonal <- row(Sigma[[1L]]) != col(Sigma[[1L]])
  for (j in seq_along(Sigma)) {
    if (any(Sigma[[j]][off_diagonal] != 0)) {
      return(FALSE)
    }
  }
  TRUE
}

spd_inverse <- function(x) {
  # Invert a symmetric positive definite matrix. `chol()` finds a stable square
  # root-like factor of the matrix, and `chol2inv()` uses that factor to compute
  # the inverse more reliably than a direct inverse call.
  chol2inv(chol(symmetrize_matrix(x)))
}

is_scalar_logical <- function(x) {
  # TRUE only for a single non-missing TRUE/FALSE value.
  is.logical(x) && length(x) == 1L && !is.na(x)
}
