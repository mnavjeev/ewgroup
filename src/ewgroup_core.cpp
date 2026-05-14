#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

inline int idx(const int row, const int col, const int ncol) {
  return row * ncol + col;
}

inline int score_idx(const int j,
                     const int k,
                     const int col,
                     const int J,
                     const int d) {
  return (j * J + k) * d + col;
}

std::vector<double> identity_matrix(const int n) {
  std::vector<double> out(n * n, 0.0);
  for (int i = 0; i < n; ++i) {
    out[idx(i, i, n)] = 1.0;
  }
  return out;
}

std::vector<double> get_sigma_slice(const Rcpp::NumericVector& Sigma_hat,
                                    const int d,
                                    const int k) {
  std::vector<double> out(d * d);
  const int offset = k * d * d;
  for (int row = 0; row < d; ++row) {
    for (int col = 0; col < d; ++col) {
      const double value = 0.5 * (
        Sigma_hat[offset + row + col * d] +
          Sigma_hat[offset + col + row * d]
      );
      out[idx(row, col, d)] = value;
    }
  }
  return out;
}

std::vector<double> invert_matrix(std::vector<double> A,
                                  const int n,
                                  const char* message) {
  std::vector<double> inv = identity_matrix(n);

  for (int col = 0; col < n; ++col) {
    int pivot = col;
    double pivot_abs = std::abs(A[idx(col, col, n)]);
    for (int row = col + 1; row < n; ++row) {
      const double candidate = std::abs(A[idx(row, col, n)]);
      if (candidate > pivot_abs) {
        pivot = row;
        pivot_abs = candidate;
      }
    }
    if (pivot_abs <= std::numeric_limits<double>::epsilon()) {
      Rcpp::stop(message);
    }
    if (pivot != col) {
      for (int j = 0; j < n; ++j) {
        std::swap(A[idx(col, j, n)], A[idx(pivot, j, n)]);
        std::swap(inv[idx(col, j, n)], inv[idx(pivot, j, n)]);
      }
    }

    const double scale = A[idx(col, col, n)];
    for (int j = 0; j < n; ++j) {
      A[idx(col, j, n)] /= scale;
      inv[idx(col, j, n)] /= scale;
    }

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

std::vector<double> inverse_shifted_covariance(const std::vector<double>& S,
                                               const int d,
                                               const double gamma,
                                               const char* message) {
  std::vector<double> M = identity_matrix(d);
  for (int i = 0; i < d * d; ++i) {
    M[i] -= gamma * S[i];
  }
  return invert_matrix(M, d, message);
}

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

std::vector<double> mat_mat(const std::vector<double>& A,
                            const std::vector<double>& B,
                            const int d) {
  std::vector<double> out(d * d, 0.0);
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

double matrix_trace(const std::vector<double>& A, const int d) {
  double out = 0.0;
  for (int i = 0; i < d; ++i) {
    out += A[idx(i, i, d)];
  }
  return out;
}

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

// [[Rcpp::export]]
Rcpp::List ewgroup_core_cpp(const Rcpp::NumericMatrix& B,
                            const Rcpp::NumericVector& Sigma_hat,
                            const double sigma2,
                            const double gamma,
                            const bool return_jacobian = true,
                            const bool return_weights = true) {
  Rcpp::IntegerVector dims = Sigma_hat.attr("dim");
  if (dims.size() != 3) {
    Rcpp::stop("Sigma_hat must be a d by d by J array.");
  }

  const int J = B.nrow();
  const int d = B.ncol();
  if (dims[0] != d || dims[1] != d || dims[2] != J) {
    Rcpp::stop("Sigma_hat dimensions must match beta_hat.");
  }

  std::vector< std::vector<double> > Sigma(J);
  std::vector< std::vector<double> > Omega(J);
  std::vector<double> trace_sigma(J);

  for (int k = 0; k < J; ++k) {
    Sigma[k] = get_sigma_slice(Sigma_hat, d, k);
    trace_sigma[k] = matrix_trace(Sigma[k], d);
    Omega[k] = inverse_shifted_covariance(
      Sigma[k],
      d,
      gamma,
      "gamma times the largest covariance eigenvalue must be below one."
    );
    for (int i = 0; i < d * d; ++i) {
      Omega[k][i] *= gamma;
    }
  }

  if (!return_jacobian && !return_weights) {
    Rcpp::NumericMatrix tilde(J, d);
    std::vector<double> trace_terms(J, 0.0);
    std::vector<double> derivative_trace(J, 0.0);
    std::vector<double> log_weights(J, 0.0);
    std::vector<double> weights(J, 0.0);
    std::vector<double> score_row(J * d, 0.0);
    std::vector<double> diff(d);
    std::vector<double> transformed(d);

    for (int j = 0; j < J; ++j) {
      double row_max = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < J; ++k) {
        for (int col = 0; col < d; ++col) {
          diff[col] = B(k, col) - B(j, col);
        }
        mat_vec_into(Omega[k], diff, d, transformed);
        double quad = 0.0;
        for (int col = 0; col < d; ++col) {
          quad += diff[col] * transformed[col];
          score_row[idx(k, col, d)] = transformed[col] / sigma2;
        }
        log_weights[k] = -0.5 * quad / sigma2;
        row_max = std::max(row_max, log_weights[k]);
      }

      double row_sum = 0.0;
      for (int k = 0; k < J; ++k) {
        weights[k] = std::exp(log_weights[k] - row_max);
        row_sum += weights[k];
      }

      std::vector<double> Gamma(d * d, 0.0);
      std::vector<double> weighted_mean(d, 0.0);
      std::vector<double> score_bar(d, 0.0);
      for (int k = 0; k < J; ++k) {
        weights[k] /= row_sum;
        const double w = weights[k];
        for (int row = 0; row < d; ++row) {
          weighted_mean[row] += w * B(k, row);
          score_bar[row] += w * score_row[idx(k, row, d)];
          for (int col = 0; col < d; ++col) {
            Gamma[idx(row, col, d)] += gamma * w * Sigma[k][idx(row, col, d)];
          }
        }
      }

      std::vector<double> H = inverse_shifted_covariance(
        Gamma,
        d,
        1.0,
        "Encountered non-invertible debiasing matrix."
      );
      std::vector<double> N(d);
      for (int col = 0; col < d; ++col) {
        N[col] = weighted_mean[col] - B(j, col);
      }
      std::vector<double> s = mat_vec(H, N, d);
      for (int col = 0; col < d; ++col) {
        tilde(j, col) = B(j, col) + s[col];
      }

      std::vector<double> Sigma_s_all(J * d, 0.0);
      for (int k = 0; k < J; ++k) {
        mat_vec_into(Sigma[k], s, d, transformed);
        for (int col = 0; col < d; ++col) {
          Sigma_s_all[idx(k, col, d)] = transformed[col];
        }
      }

      std::vector<double> Sigma_H = mat_mat(Sigma[j], H, d);
      double trace_update =
        -(1.0 - weights[j]) * matrix_trace(Sigma_H, d);
      for (int k = 0; k < J; ++k) {
        for (int col = 0; col < d; ++col) {
          diff[col] = B(k, col) - B(j, col) +
            gamma * Sigma_s_all[idx(k, col, d)];
        }
        mat_vec_into(Sigma_H, diff, d, transformed);
        for (int col = 0; col < d; ++col) {
          const double q = weights[k] * (score_row[idx(k, col, d)] - score_bar[col]);
          trace_update += q * transformed[col];
        }
      }

      trace_terms[j] = trace_update;
      derivative_trace[j] = trace_sigma[j] + trace_terms[j];
    }

    double sure_A = 0.0;
    for (int j = 0; j < J; ++j) {
      for (int col = 0; col < d; ++col) {
        const double residual = tilde(j, col) - B(j, col);
        sure_A += residual * residual;
      }
    }

    double sure_D_sum = 0.0;
    for (int j = 0; j < J; ++j) {
      sure_D_sum += trace_terms[j];
    }
    const double sure_D = sigma2 * sure_D_sum;
    const double alpha_unconstrained =
      sure_A <= std::numeric_limits<double>::epsilon() ? 0.0 : -sure_D / sure_A;
    const double alpha = std::min(1.0, std::max(0.0, alpha_unconstrained));

    return Rcpp::List::create(
      Rcpp::Named("tilde") = tilde,
      Rcpp::Named("alpha") = alpha,
      Rcpp::Named("alpha_unconstrained") = alpha_unconstrained,
      Rcpp::Named("sure_A") = sure_A,
      Rcpp::Named("sure_D") = sure_D,
      Rcpp::Named("weights") = R_NilValue,
      Rcpp::Named("jacobian") = Rcpp::List(0),
      Rcpp::Named("derivative_trace") = Rcpp::wrap(derivative_trace)
    );
  }

  std::vector<double> weights(J * J, 0.0);
  std::vector<double> log_weights(J * J, 0.0);
  std::vector<double> scores(J * J * d, 0.0);
  std::vector<double> diff(d);
  std::vector<double> transformed(d);

  for (int k = 0; k < J; ++k) {
    for (int j = 0; j < J; ++j) {
      for (int col = 0; col < d; ++col) {
        diff[col] = B(k, col) - B(j, col);
      }
      mat_vec_into(Omega[k], diff, d, transformed);
      double quad = 0.0;
      for (int col = 0; col < d; ++col) {
        quad += diff[col] * transformed[col];
        scores[score_idx(j, k, col, J, d)] = transformed[col] / sigma2;
      }
      log_weights[idx(j, k, J)] = -0.5 * quad / sigma2;
    }
  }

  for (int j = 0; j < J; ++j) {
    double row_max = log_weights[idx(j, 0, J)];
    for (int k = 1; k < J; ++k) {
      row_max = std::max(row_max, log_weights[idx(j, k, J)]);
    }
    double row_sum = 0.0;
    for (int k = 0; k < J; ++k) {
      const double value = std::exp(log_weights[idx(j, k, J)] - row_max);
      weights[idx(j, k, J)] = value;
      row_sum += value;
    }
    for (int k = 0; k < J; ++k) {
      weights[idx(j, k, J)] /= row_sum;
    }
  }

  Rcpp::NumericMatrix tilde(J, d);
  Rcpp::NumericMatrix weights_out;
  if (return_weights) {
    weights_out = Rcpp::NumericMatrix(J, J);
  }
  std::vector<double> trace_terms(J, 0.0);
  std::vector<double> derivative_trace(J, 0.0);
  Rcpp::List jacobian(return_jacobian ? J : 0);

  for (int j = 0; j < J; ++j) {
    std::vector<double> Gamma(d * d, 0.0);
    std::vector<double> weighted_mean(d, 0.0);
    for (int k = 0; k < J; ++k) {
      const double w = weights[idx(j, k, J)];
      if (return_weights) {
        weights_out(j, k) = w;
      }
      for (int row = 0; row < d; ++row) {
        weighted_mean[row] += w * B(k, row);
        for (int col = 0; col < d; ++col) {
          Gamma[idx(row, col, d)] += gamma * w * Sigma[k][idx(row, col, d)];
        }
      }
    }

    std::vector<double> H = inverse_shifted_covariance(
      Gamma,
      d,
      1.0,
      "Encountered non-invertible debiasing matrix."
    );
    std::vector<double> N(d);
    for (int col = 0; col < d; ++col) {
      N[col] = weighted_mean[col] - B(j, col);
    }
    std::vector<double> s = mat_vec(H, N, d);
    for (int col = 0; col < d; ++col) {
      tilde(j, col) = B(j, col) + s[col];
    }

    std::vector<double> score_bar(d, 0.0);
    for (int k = 0; k < J; ++k) {
      const double w = weights[idx(j, k, J)];
      for (int col = 0; col < d; ++col) {
        score_bar[col] += w * scores[score_idx(j, k, col, J, d)];
      }
    }

    std::vector<double> Sigma_s_all(J * d, 0.0);
    for (int k = 0; k < J; ++k) {
      mat_vec_into(Sigma[k], s, d, transformed);
      for (int col = 0; col < d; ++col) {
        Sigma_s_all[idx(k, col, d)] = transformed[col];
      }
    }

    if (!return_jacobian) {
      std::vector<double> Sigma_H = mat_mat(Sigma[j], H, d);
      // Trace-only form of tr(Sigma_j (P_j - I)) avoids materializing P_j.
      double trace_update =
        -(1.0 - weights[idx(j, j, J)]) * matrix_trace(Sigma_H, d);
      for (int k = 0; k < J; ++k) {
        for (int col = 0; col < d; ++col) {
          diff[col] = B(k, col) - B(j, col) +
            gamma * Sigma_s_all[idx(k, col, d)];
        }
        mat_vec_into(Sigma_H, diff, d, transformed);
        for (int col = 0; col < d; ++col) {
          const double q =
            weights[idx(j, k, J)] *
            (scores[score_idx(j, k, col, J, d)] - score_bar[col]);
          trace_update += q * transformed[col];
        }
      }

      trace_terms[j] = trace_update;
      derivative_trace[j] = trace_sigma[j] + trace_terms[j];
      continue;
    }

    Rcpp::NumericMatrix P(d, d);
    std::vector<double> P_vec(d * d, 0.0);
    for (int r = 0; r < d; ++r) {
      std::vector<double> dN(d, 0.0);
      std::vector<double> dGamma_s(d, 0.0);

      for (int k = 0; k < J; ++k) {
        const double dw =
          weights[idx(j, k, J)] *
          (scores[score_idx(j, k, r, J, d)] - score_bar[r]);
        for (int col = 0; col < d; ++col) {
          dN[col] += dw * (B(k, col) - B(j, col));
        }
        for (int row = 0; row < d; ++row) {
          dGamma_s[row] += gamma * dw * Sigma_s_all[idx(k, row, d)];
        }
      }
      dN[r] -= 1.0 - weights[idx(j, j, J)];

      std::vector<double> H_dGamma_s = mat_vec(H, dGamma_s, d);
      std::vector<double> H_dN = mat_vec(H, dN, d);
      for (int row = 0; row < d; ++row) {
        double value = H_dGamma_s[row] + H_dN[row];
        if (row == r) {
          value += 1.0;
        }
        P(row, r) = value;
        P_vec[idx(row, r, d)] = value;
      }
    }

    std::vector<double> P_minus_I = P_vec;
    for (int i = 0; i < d; ++i) {
      P_minus_I[idx(i, i, d)] -= 1.0;
    }
    trace_terms[j] = trace_product(Sigma[j], P_minus_I, d);
    derivative_trace[j] = trace_sigma[j] + trace_terms[j];
    if (return_jacobian) {
      jacobian[j] = P;
    }
  }

  double sure_A = 0.0;
  for (int j = 0; j < J; ++j) {
    for (int col = 0; col < d; ++col) {
      const double residual = tilde(j, col) - B(j, col);
      sure_A += residual * residual;
    }
  }

  double sure_D_sum = 0.0;
  for (int j = 0; j < J; ++j) {
    sure_D_sum += trace_terms[j];
  }
  const double sure_D = sigma2 * sure_D_sum;
  const double alpha_unconstrained =
    sure_A <= std::numeric_limits<double>::epsilon() ? 0.0 : -sure_D / sure_A;
  const double alpha = std::min(1.0, std::max(0.0, alpha_unconstrained));
  Rcpp::RObject weights_return = R_NilValue;
  if (return_weights) {
    weights_return = weights_out;
  }

  return Rcpp::List::create(
    Rcpp::Named("tilde") = tilde,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_unconstrained") = alpha_unconstrained,
    Rcpp::Named("sure_A") = sure_A,
    Rcpp::Named("sure_D") = sure_D,
    Rcpp::Named("weights") = weights_return,
    Rcpp::Named("jacobian") = jacobian,
    Rcpp::Named("derivative_trace") = Rcpp::wrap(derivative_trace)
  );
}
