# A deliberately direct implementation of the paper's matrix formulas for
# small inputs. Its dense weights and finite differences are independent of
# the streaming analytic derivative used by the package.
paper_tilde <- function(B, S, sigma2, gamma) {
  J <- nrow(B)
  d <- ncol(B)
  weights <- matrix(0, J, J)
  result <- B
  for (j in seq_len(J)) {
    for (k in seq_len(J)) {
      z <- B[j, ] - B[k, ]
      weights[j, k] <- exp(-gamma * sum(z * solve(diag(d) - gamma * S[[k]], z)) /
                            (2 * sigma2))
    }
    weights[j, ] <- weights[j, ] / sum(weights[j, ])
    G <- Reduce(`+`, Map(function(w, x) w * gamma * x, weights[j, ], S))
    result[j, ] <- solve(diag(d) - G,
                        colSums(weights[j, ] * B) - G %*% B[j, ])
  }
  list(tilde = result, weights = weights)
}

repair_inputs <- function(kind) {
  if (kind == "scalar") {
    return(list(B = matrix(c(-1.1, -0.95, 0.2, 0.27, 1.3)),
                S = lapply(c(0.8, 1.1, 0.9, 1.2, 1), function(x) matrix(x))))
  }
  B <- rbind(c(-1, 0.25), c(-0.75, 0), c(0.5, -0.25), c(0.75, -0.5))
  S <- lapply(seq_len(nrow(B)), function(k) {
    x <- diag(c(0.7 + k / 10, 1.4 - k / 20))
    if (kind == "full") x[1, 2] <- x[2, 1] <- 0.1 + k / 100
    x
  })
  list(B = B, S = S)
}

test_that("compiled and fallback modes match paper formulas and numerical derivatives", {
  for (kind in c("scalar", "diagonal", "full")) {
    x <- repair_inputs(kind)
    B <- x$B
    S <- x$S
    sigma2 <- 0.1
    gamma <- 0.1
    ref <- paper_tilde(B, S, sigma2, gamma)
    derivatives <- lapply(seq_len(nrow(B)), function(j) {
      P <- matrix(0, ncol(B), ncol(B))
      for (r in seq_len(ncol(B))) {
        plus <- minus <- B
        plus[j, r] <- plus[j, r] + 1e-5
        minus[j, r] <- minus[j, r] - 1e-5
        P[, r] <- (paper_tilde(plus, S, sigma2, gamma)$tilde[j, ] -
                     paper_tilde(minus, S, sigma2, gamma)$tilde[j, ]) / 2e-5
      }
      P
    })
    D <- sigma2 * sum(vapply(seq_len(nrow(B)), function(j) {
      sum(diag(S[[j]] %*% (derivatives[[j]] - diag(ncol(B)))))
    }, numeric(1)))
    ref_alpha <- max(0, min(1, -D / sum((ref$tilde - B)^2)))
    for (derivative in c(FALSE, TRUE)) for (weights in c(FALSE, TRUE)) {
      cpp <- ewgroup:::ewgroup_core(B, S, sigma2, gamma, derivative, weights)
      pure <- ewgroup:::ewgroup_core_r(B, S, sigma2, gamma, derivative, weights)
      expect_equal(cpp, pure, tolerance = 2e-12)
      expect_equal(cpp$tilde, ref$tilde, tolerance = 1e-12)
      expect_equal(cpp$alpha, ref_alpha, tolerance = 2e-8)
      if (weights) expect_equal(cpp$weights, ref$weights, tolerance = 1e-12)
      else expect_null(pure$weights)
      if (derivative) expect_equal(cpp$jacobian, derivatives, tolerance = 2e-8)
      else expect_length(pure$jacobian, 0)
    }
  }
})

test_that("unit rescaling preserves the SURE decision and estimator", {
  for (kind in c("scalar", "diagonal", "full")) {
    x <- repair_inputs(kind)
    for (core in list(ewgroup:::ewgroup_core, ewgroup:::ewgroup_core_r)) {
      base <- core(x$B, x$S, 0.1, 0.1, TRUE, TRUE)
      for (scale in c(1e-110, 1e-8, 1e8, 1e110)) {
        scaled <- core(x$B * scale, lapply(x$S, function(S) S * scale^2),
                       0.1, 0.1 / scale^2, TRUE, TRUE)
        expect_equal(scaled$alpha, base$alpha, tolerance = 2e-12)
        expect_equal(scaled$tilde / scale, base$tilde, tolerance = 2e-12)
        expect_equal(scaled$weights, base$weights, tolerance = 2e-12)
        expect_equal(scaled$jacobian, base$jacobian, tolerance = 2e-12)
        expect_equal(scaled$sure_A / scale^2, base$sure_A, tolerance = 2e-12)
        expect_equal(scaled$sure_D / scale^2, base$sure_D, tolerance = 2e-12)
      }
    }
    base <- ewgroup(x$B, x$S, 0.1)
    small <- ewgroup(x$B * 1e-8, lapply(x$S, function(S) S * 1e-16), 0.1)
    expect_equal(small$alpha, base$alpha, tolerance = 2e-12)
  }
})

test_that("common shifts and cell ordering preserve smoothing decisions", {
  x <- repair_inputs("full")
  for (core in list(ewgroup:::ewgroup_core, ewgroup:::ewgroup_core_r)) {
    fit <- core(x$B, x$S, 0.1, 0.1, TRUE, TRUE)
    shifted <- core(x$B + 2^40, x$S, 0.1, 0.1, TRUE, TRUE)
    expect_equal(shifted$adjustment, fit$adjustment, tolerance = 1e-12)
    expect_equal(shifted$alpha, fit$alpha, tolerance = 1e-12)
    expect_equal(shifted$sure_A, fit$sure_A, tolerance = 1e-12)
    order <- c(3, 1, 4, 2)
    reordered <- core(x$B[order, ], x$S[order], 0.1, 0.1, TRUE, TRUE)
    expect_equal(reordered$tilde, fit$tilde[order, ], tolerance = 1e-12)
    expect_equal(reordered$weights, fit$weights[order, order], tolerance = 1e-12)
    expect_equal(reordered$alpha, fit$alpha, tolerance = 1e-12)
  }
})

test_that("zero displacement has alpha zero and small displacement remains usable", {
  for (core in list(ewgroup:::ewgroup_core, ewgroup:::ewgroup_core_r)) {
    for (J in c(1L, 4L)) {
      fit <- core(matrix(rep(1, J)), rep(list(matrix(1)), J), 0.1, 0.2, TRUE, TRUE)
      expect_equal(fit$alpha, 0)
      expect_equal(fit$sure_A, 0)
      expect_equal(fit$tilde, matrix(rep(1, J)))
    }
  }
  # The squared norm underflows, but its logarithm still gives the right ratio.
  mix <- ewgroup:::sure_mixing(matrix(1e-200), -1e-200, 1e-200)
  expect_equal(mix$sure_A, 0)
  expect_equal(mix$sure_D, 0)
  expect_equal(mix$alpha, 1, tolerance = 1e-12)
})

test_that("covariance validation is nonempty, dimensionally consistent and relative", {
  expect_error(ewgroup(numeric(), numeric(), 1), "at least one")
  for (empty in list(numeric(), list(), list(matrix(numeric(), 0, 0)),
                     array(numeric(), c(2, 2, 0)), array(numeric(), c(0, 0, 1)))) {
    expect_error(ewgroup_gamma(empty), "at least one|nonempty")
  }
  expect_error(ewgroup_gamma(list(diag(2), diag(3))), "dimensions")
  for (scale in c(1e-16, 1, 1e16)) {
    asym <- matrix(c(1, 0, 0.1, 1), 2) * scale
    indefinite <- matrix(c(1, 2, 2, 1), 2) * scale
    expect_error(ewgroup_gamma(list(asym)), "symmetric")
    expect_error(ewgroup_gamma(list(indefinite)), "positive semidefinite")
    expect_error(ewgroup_gamma(c(scale, -scale * 1e-12)), "positive semidefinite")
    expect_error(ewgroup_gamma(list(diag(c(1, -1e-12)) * scale)), "positive semidefinite")
    tiny_negative <- matrix(c(1, 1 + 1e-10, 1 + 1e-10, 1), 2) * scale
    checked <- ewgroup:::normalize_covariances(list(tiny_negative))$matrices[[1]]
    expect_gte(min(eigen(checked / scale, symmetric = TRUE)$values), -1e-14)
  }
  expect_error(ewgroup_gamma(c(0, 0)), "positive eigenvalue")
  expect_equal(ewgroup(c(1, 2), c(0, 0), 1, gamma = 0.1)$alpha, 0)
})

test_that("tiny nonzero covariance entries are retained by fallback dispatch", {
  S <- rep(list(matrix(c(1, 0.25, 0.25, 1), 2) * 1e-16), 3)
  expect_false(ewgroup:::covariance_list_is_diagonal(S))
  B <- rbind(c(0, 1), c(0.2, 0.9), c(1, 0)) * 1e-8
  expect_equal(ewgroup:::ewgroup_core(B, S, 0.1, 1e15),
               ewgroup:::ewgroup_core_r(B, S, 0.1, 1e15), tolerance = 2e-12)
})

test_that("fallback dispatch passes diagnostic flags through", {
  # Source the implementation without its compiled wrapper, as fixture scripts
  # and development sessions do. No installed C++ symbol exists in this scope.
  env <- new.env(parent = baseenv())
  env$coef <- stats::coef
  env$fitted <- stats::fitted
  env$eigen <- base::eigen
  environment_core <- environment(ewgroup:::ewgroup_core)
  for (name in ls(environment_core, all.names = TRUE)) {
    if (name == "ewgroup_core_cpp") next
    value <- get(name, environment_core)
    if (is.function(value) && !is.primitive(value)) environment(value) <- env
    assign(name, value, env)
  }
  fit <- env$ewgroup(c(0, 1), c(1, 1), sigma2 = 0.1)
  expect_null(fit$weights)
  expect_null(fit$derivative)
  core <- env$ewgroup_core(matrix(c(0, 1)), list(matrix(1), matrix(1)),
                          0.1, 0.1, FALSE, FALSE)
  expect_null(core$weights)
  expect_length(core$jacobian, 0)
})

test_that("nonfinite numerical intermediates fail with actionable errors", {
  expect_error(ewgroup(c(-1e308, 1e308), c(1, 1), 1), "Nonfinite numerical")
  expect_error(ewgroup(c(0, 1e200), c(1, 1), 1), "Nonfinite numerical")
})

test_that("SURE exercises interior and both constrained endpoints", {
  cases <- list(
    list(b = c(-1.2, .4, -.3, -.5, 1, -.2), v = c(1.6, 1.9, .6, 1.8, .9, 1.2)),
    list(b = c(2.1, -1.5, .3, 1.5, -.4, .6), v = c(1.25, 1.57, 1.69, 1.84, .26, 1.65)),
    list(b = c(0, .01), v = c(1, 1))
  )
  for (i in seq_along(cases)) {
    x <- cases[[i]]
    S <- lapply(x$v, matrix)
    for (core in list(ewgroup:::ewgroup_core, ewgroup:::ewgroup_core_r)) {
      fit <- core(matrix(x$b), S, 0.1, 0.1, TRUE, FALSE)
      expect_gt(fit$sure_A, 0)
      if (i == 1L) {
        expect_gt(fit$alpha, 0)
        expect_lt(fit$alpha, 1)
      } else expect_equal(fit$alpha, if (i == 2L) 0 else 1)
      scaled <- core(matrix(x$b * 1e-8), lapply(S, function(s) s * 1e-16),
                     0.1, 1e15, TRUE, FALSE)
      expect_equal(scaled$alpha, fit$alpha, tolerance = 2e-12)
    }
  }
})

test_that("fallback without weights avoids quadratic allocations", {
  skip_if_not(capabilities("profmem"))
  J <- 300L
  B <- matrix(seq(-1, 1, length.out = J))
  S <- rep(list(matrix(1)), J)
  # Warm up lazy loads before recording allocations made by the estimator.
  invisible(ewgroup:::ewgroup_core_r(B[1:2, , drop = FALSE], S[1:2], .1, .1, TRUE, FALSE))
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  Rprofmem(path)
  tryCatch(ewgroup:::ewgroup_core_r(B, S, .1, .1, TRUE, FALSE),
           finally = Rprofmem(NULL))
  lines <- readLines(path)
  allocations <- as.numeric(sub(" .*", "", lines[grepl("^[0-9]+", lines)]))
  expect_lt(max(allocations), 8 * J^2)
})
