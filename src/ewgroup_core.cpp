#include <Rcpp.h>

// Rcpp is the bridge between R and C++. It supplies C++ versions of common R
// objects, such as NumericMatrix and List, and it lets R call the function at
// the bottom of this file through the generated wrapper files.
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

// Store a row-and-column location as one position in a flat vector. The code
// keeps many small matrices in plain vectors because that is faster in C++.
inline size_t idx(const int row, const int col, const int ncol) {
  return static_cast<size_t>(row) * ncol + col;
}

// Check products before allocating or indexing. R arrays and C++ vectors must
// both be able to represent the requested storage size.
size_t checked_size(const int rows, const int cols) {
  const size_t a = static_cast<size_t>(rows);
  const size_t b = static_cast<size_t>(cols);
  if (rows < 1 || cols < 1 || a > std::numeric_limits<size_t>::max() / b ||
      a * b > std::vector<double>().max_size() ||
      a * b > static_cast<size_t>(R_XLEN_T_MAX)) {
    Rcpp::stop("Requested estimator workspace is too large.");
  }
  return a * b;
}

void require_finite(const double x) {
  if (!std::isfinite(x)) {
    Rcpp::stop("Nonfinite numerical result; rescale beta_hat and its covariance.");
  }
}

// Make an identity matrix: ones on the diagonal and zeros elsewhere. Multiplying
// by an identity matrix leaves a vector unchanged, so it plays the role of "1"
// for matrix formulas.
std::vector<double> identity_matrix(const int n) {
  std::vector<double> out(checked_size(n, n), 0.0);
  for (int i = 0; i < n; ++i) {
    out[idx(i, i, n)] = 1.0;
  }
  return out;
}

// Pull one covariance matrix out of R's d by d by J array. R stores arrays in
// column order, so this helper reads the R layout and writes the C++ row order.
// It also averages mirrored entries to remove tiny roundoff differences.
std::vector<double> get_sigma_slice(const Rcpp::NumericVector& Sigma_hat,
                                    const int d,
                                    const int k) {
  std::vector<double> out(checked_size(d, d));
  const size_t offset = static_cast<size_t>(k) * d * d;
  for (int row = 0; row < d; ++row) {
    for (int col = 0; col < d; ++col) {
      const double value =
        0.5 * Sigma_hat[offset + row + static_cast<size_t>(col) * d] +
        0.5 * Sigma_hat[offset + col + static_cast<size_t>(row) * d];
      out[idx(row, col, d)] = value;
    }
  }
  return out;
}

// Invert a small square matrix using row operations. Inverting means finding a
// matrix that undoes the original matrix. The `message` text is shown if the
// matrix cannot be safely inverted.
std::vector<double> invert_matrix(std::vector<double> A,
                                  const int n,
                                  const char* message) {
  // Start with the identity matrix. The same row operations applied to `A` are
  // also applied to `inv`; when `A` becomes the identity matrix, `inv` has
  // become the inverse.
  std::vector<double> inv = identity_matrix(n);

  for (int col = 0; col < n; ++col) {
    // Choose the row with the largest available entry in this column. This
    // reduces avoidable numerical problems when dividing below.
    int pivot = col;
    double pivot_abs = std::abs(A[idx(col, col, n)]);
    for (int row = col + 1; row < n; ++row) {
      const double candidate = std::abs(A[idx(row, col, n)]);
      if (candidate > pivot_abs) {
        pivot = row;
        pivot_abs = candidate;
      }
    }
    if (!std::isfinite(pivot_abs) || pivot_abs == 0.0) {
      Rcpp::stop(message);
    }

    // Move the best pivot row into the current row.
    if (pivot != col) {
      for (int j = 0; j < n; ++j) {
        std::swap(A[idx(col, j, n)], A[idx(pivot, j, n)]);
        std::swap(inv[idx(col, j, n)], inv[idx(pivot, j, n)]);
      }
    }

    // Scale the pivot row so the pivot entry equals one.
    const double scale = A[idx(col, col, n)];
    for (int j = 0; j < n; ++j) {
      A[idx(col, j, n)] /= scale;
      inv[idx(col, j, n)] /= scale;
    }

    // Use the pivot row to clear this column from every other row.
    for (int row = 0; row < n; ++row) {
      if (row == col) {
        continue;
      }
      const double factor = A[idx(row, col, n)];
      if (factor == 0.0) {
        continue;
      }
      for (int j = 0; j < n; ++j) {
        A[idx(row, j, n)] -= factor * A[idx(col, j, n)];
        inv[idx(row, j, n)] -= factor * inv[idx(col, j, n)];
      }
    }
  }

  return inv;
}

// Build and invert I - gamma * S, where S is a covariance matrix. This matrix
// appears repeatedly in the estimator's weighting and correction formulas.
std::vector<double> inverse_shifted_covariance(const std::vector<double>& S,
                                               const int d,
                                               const double gamma,
                                               const char* message) {
  std::vector<double> M = identity_matrix(d);
  for (size_t i = 0; i < checked_size(d, d); ++i) {
    M[i] -= gamma * S[i];
  }
  return invert_matrix(M, d, message);
}

// Multiply a matrix by a vector and return the resulting vector.
std::vector<double> mat_vec(const std::vector<double>& A,
                            const std::vector<double>& x,
                            const int d) {
  std::vector<double> out(d, 0.0);
  for (int row = 0; row < d; ++row) {
    double value = 0.0;
    for (int col = 0; col < d; ++col) {
      value += A[idx(row, col, d)] * x[col];
    }
    out[row] = value;
  }
  return out;
}

// Same matrix-times-vector calculation as `mat_vec()`, but writes into an
// existing output vector. Reusing the output vector avoids repeated memory
// allocation inside large loops.
void mat_vec_into(const std::vector<double>& A,
                  const std::vector<double>& x,
                  const int d,
                  std::vector<double>& out) {
  for (int row = 0; row < d; ++row) {
    double value = 0.0;
    for (int col = 0; col < d; ++col) {
      value += A[idx(row, col, d)] * x[col];
    }
    out[row] = value;
  }
}

// Multiply two square matrices.
std::vector<double> mat_mat(const std::vector<double>& A,
                            const std::vector<double>& B,
                            const int d) {
  std::vector<double> out(checked_size(d, d), 0.0);
  for (int row = 0; row < d; ++row) {
    for (int col = 0; col < d; ++col) {
      double value = 0.0;
      for (int inner = 0; inner < d; ++inner) {
        value += A[idx(row, inner, d)] * B[idx(inner, col, d)];
      }
      out[idx(row, col, d)] = value;
    }
  }
  return out;
}

// Add up the diagonal entries of a matrix.
double matrix_trace(const std::vector<double>& A, const int d) {
  double out = 0.0;
  for (int i = 0; i < d; ++i) {
    out += A[idx(i, i, d)];
  }
  return out;
}

// Compute the trace of A times B without first building the product matrix.
// This is faster and uses less memory for the SURE correction.
double trace_product(const std::vector<double>& A,
                     const std::vector<double>& B,
                     const int d) {
  double out = 0.0;
  for (int i = 0; i < d; ++i) {
    for (int ell = 0; ell < d; ++ell) {
      out += A[idx(i, ell, d)] * B[idx(ell, i, d)];
    }
  }
  return out;
}

} // namespace

// The Rcpp export tag tells Rcpp to make an R-callable wrapper for this C++
// function. Users do not call this function directly; ewgroup() calls it after
// checking and reshaping the inputs.
// [[Rcpp::export]]
Rcpp::List ewgroup_core_cpp(const Rcpp::NumericMatrix& B,
                            const Rcpp::NumericVector& Sigma_hat,
                            const double sigma2,
                            const double gamma,
                            const bool return_jacobian = true,
                            const bool return_weights = true) {
  // `Sigma_hat` arrives from R as a numeric vector with array dimensions
  // attached. Confirm that it is a three-dimensional covariance array.
  Rcpp::IntegerVector dims = Sigma_hat.attr("dim");
  if (dims.size() != 3) {
    Rcpp::stop("Sigma_hat must be a d by d by J array.");
  }

  // J is the number of cells. d is the number of coefficients per cell.
  const int J = B.nrow();
  const int d = B.ncol();
  if (dims[0] != d || dims[1] != d || dims[2] != J) {
    Rcpp::stop("Sigma_hat dimensions must match beta_hat.");
  }

  checked_size(J, d);
  const size_t dd = checked_size(d, d);
  if (static_cast<size_t>(J) > static_cast<size_t>(R_XLEN_T_MAX) / dd) {
    Rcpp::stop("Requested covariance workspace is too large.");
  }
  if (!std::isfinite(sigma2) || sigma2 <= 0 ||
      !std::isfinite(gamma) || gamma <= 0) {
    Rcpp::stop("sigma2 and gamma must be positive finite scalars.");
  }
  std::vector< std::vector<double> > Sigma(J), scaled_sigma(J), Omega(J);
  std::vector<double> trace_sigma(J);
  bool diagonal = true;
  for (int k = 0; k < J; ++k) {
    if (k % 256 == 0) Rcpp::checkUserInterrupt();
    Sigma[k] = get_sigma_slice(Sigma_hat, d, k);
    scaled_sigma[k].resize(dd);
    for (int row = 0; row < d; ++row) {
      require_finite(B(k, row));
      for (int col = 0; col < d; ++col) {
        const size_t i = idx(row, col, d);
        require_finite(Sigma[k][i]);
        scaled_sigma[k][i] = gamma * Sigma[k][i];
        require_finite(scaled_sigma[k][i]);
        if (row != col && Sigma[k][i] != 0.0) diagonal = false;
      }
    }
    trace_sigma[k] = matrix_trace(Sigma[k], d);
  }
  for (int k = 0; k < J; ++k) {
    if (diagonal) {
      Omega[k].assign(dd, 0.0);
      for (int r = 0; r < d; ++r) {
        const size_t i = idx(r, r, d);
        const double denom = 1.0 - scaled_sigma[k][i];
        if (denom <= 0) Rcpp::stop("Non-positive distance-scaling matrix.");
        Omega[k][i] = gamma / denom;
        require_finite(Omega[k][i]);
      }
    } else {
      Omega[k] = inverse_shifted_covariance(
        scaled_sigma[k], d, 1.0, "Non-invertible distance-scaling matrix.");
      for (size_t i = 0; i < dd; ++i) {
        Omega[k][i] *= gamma;
        require_finite(Omega[k][i]);
      }
    }
  }

  // All modes stream one target cell. Only the explicitly requested weight
  // matrix occupies quadratic space; Jacobians require J*d*d entries.
  Rcpp::NumericMatrix tilde(J, d), adjustment(J, d);
  Rcpp::NumericMatrix weights_out;
  if (return_weights) {
    checked_size(J, J);
    weights_out = Rcpp::NumericMatrix(J, J);
  }
  Rcpp::List jacobian(return_jacobian ? J : 0);
  std::vector<double> derivative_trace(J), trace_terms(J);
  std::vector<double> weights(J), scores(checked_size(J, d));
  std::vector<double> diff(d), transformed(d), sigma_s(d);
  double adjustment_scale = 0.0;

  for (int j = 0; j < J; ++j) {
    Rcpp::checkUserInterrupt();
    double row_max = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < J; ++k) {
      if (k % 16384 == 0) Rcpp::checkUserInterrupt();
      for (int r = 0; r < d; ++r) {
        diff[r] = B(k, r) - B(j, r);
        require_finite(diff[r]);
      }
      if (diagonal) {
        for (int r = 0; r < d; ++r) transformed[r] = Omega[k][idx(r, r, d)] * diff[r];
      } else {
        mat_vec_into(Omega[k], diff, d, transformed);
      }
      double quad = 0.0;
      for (int r = 0; r < d; ++r) {
        const double score = transformed[r] / sigma2;
        require_finite(score);
        scores[idx(k, r, d)] = score;
        quad += diff[r] * score;
      }
      require_finite(quad);
      weights[k] = -0.5 * quad;
      row_max = std::max(row_max, weights[k]);
    }
    double row_sum = 0.0;
    for (int k = 0; k < J; ++k) {
      weights[k] = std::exp(weights[k] - row_max);
      row_sum += weights[k];
    }
    require_finite(row_sum);
    std::vector<double> Gamma(dd, 0.0), N(d, 0.0), score_bar(d, 0.0);
    for (int k = 0; k < J; ++k) {
      const double w = weights[k] /= row_sum;
      if (return_weights) weights_out(j, k) = w;
      for (int r = 0; r < d; ++r) {
        // Average differences directly, so common shifts do not contaminate
        // the smoothing movement or SURE through cancellation of large means.
        N[r] += w * (B(k, r) - B(j, r));
        score_bar[r] += w * scores[idx(k, r, d)];
      }
      if (diagonal) {
        for (int r = 0; r < d; ++r) Gamma[idx(r, r, d)] += w * scaled_sigma[k][idx(r, r, d)];
      } else {
        for (size_t i = 0; i < dd; ++i) Gamma[i] += w * scaled_sigma[k][i];
      }
    }
    std::vector<double> H(dd, 0.0);
    if (diagonal) {
      for (int r = 0; r < d; ++r) {
        const size_t i = idx(r, r, d);
        const double denom = 1.0 - Gamma[i];
        if (denom <= 0) Rcpp::stop("Non-positive debiasing matrix.");
        H[i] = 1.0 / denom;
      }
    } else {
      H = inverse_shifted_covariance(Gamma, d, 1.0, "Non-invertible debiasing matrix.");
    }
    std::vector<double> s = mat_vec(H, N, d);
    for (int r = 0; r < d; ++r) {
      require_finite(s[r]);
      adjustment(j, r) = s[r];
      adjustment_scale = std::max(adjustment_scale, std::abs(s[r]));
      tilde(j, r) = B(j, r) + s[r];
      require_finite(tilde(j, r));
    }

    // U is the derivative of the displacement. Computing U directly avoids
    // subtracting an identity matrix from a nearly identical Jacobian.
    std::vector<double> U(dd, 0.0);
    for (int r = 0; r < d; ++r) U[idx(r, r, d)] = -(1.0 - weights[j]);
    for (int k = 0; k < J; ++k) {
      if (k % 16384 == 0) Rcpp::checkUserInterrupt();
      if (diagonal) {
        for (int r = 0; r < d; ++r) sigma_s[r] = scaled_sigma[k][idx(r, r, d)] * s[r];
      } else {
        mat_vec_into(scaled_sigma[k], s, d, sigma_s);
      }
      for (int r = 0; r < d; ++r) diff[r] = B(k, r) - B(j, r) + sigma_s[r];
      for (int col = 0; col < d; ++col) {
        const double dw = weights[k] * (scores[idx(k, col, d)] - score_bar[col]);
        for (int row = 0; row < d; ++row) U[idx(row, col, d)] += diff[row] * dw;
      }
    }
    U = mat_mat(H, U, d);
    for (size_t i = 0; i < dd; ++i) require_finite(U[i]);
    trace_terms[j] = trace_product(Sigma[j], U, d);
    require_finite(trace_terms[j]);
    derivative_trace[j] = trace_sigma[j] + trace_terms[j];
    if (return_jacobian) {
      Rcpp::NumericMatrix P(d, d);
      for (int row = 0; row < d; ++row) {
        for (int col = 0; col < d; ++col) P(row, col) = U[idx(row, col, d)] + (row == col);
      }
      jacobian[j] = P;
    }
  }

  // Normalize the displacement norm before forming the ratio. In particular,
  // small but nonzero adjustments are not compared to an absolute epsilon.
  double norm_sum = 0.0;
  if (adjustment_scale > 0) {
    for (int j = 0; j < J; ++j) {
      for (int r = 0; r < d; ++r) {
        const double value = adjustment(j, r) / adjustment_scale;
        norm_sum += value * value;
      }
    }
  }
  double trace_sum = 0.0;
  for (int j = 0; j < J; ++j) trace_sum += trace_terms[j];
  require_finite(trace_sum);
  const double sure_A = adjustment_scale == 0 ? 0.0 :
    std::exp(2.0 * std::log(adjustment_scale) + std::log(norm_sum));
  const double sure_D = sigma2 * trace_sum;
  require_finite(sure_A);
  require_finite(sure_D);
  double alpha_unconstrained = 0.0;
  if (adjustment_scale > 0 && trace_sum != 0) {
    const double log_ratio = std::log(sigma2) + std::log(std::abs(trace_sum)) -
      2.0 * std::log(adjustment_scale) - std::log(norm_sum);
    alpha_unconstrained = -std::copysign(std::exp(log_ratio), trace_sum);
  }
  const double alpha = std::min(1.0, std::max(0.0, alpha_unconstrained));
  Rcpp::RObject weights_return = R_NilValue;
  if (return_weights) weights_return = weights_out;
  return Rcpp::List::create(
    Rcpp::Named("tilde") = tilde,
    Rcpp::Named("adjustment") = adjustment,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_unconstrained") = alpha_unconstrained,
    Rcpp::Named("sure_A") = sure_A,
    Rcpp::Named("sure_D") = sure_D,
    Rcpp::Named("weights") = weights_return,
    Rcpp::Named("jacobian") = jacobian,
    Rcpp::Named("derivative_trace") = Rcpp::wrap(derivative_trace)
  );
}
