#include "ecdf_impl.h"

// Zero-inflated empirical distribution, one row per pixel or zone:
//   P(X = 0)  = p_zero                            (fraction of zero years)
//   P(X <= x) = p_zero + (1 - p_zero) * F(x),     x > 0
// where F is the type-7 ECDF of the positive values (R's default quantile()).
// A zero observation maps to the mid-probability p_zero / 2, as in zig_p.
// ecdf_p and ecdf_q are inverses on the positive support; p < p_zero maps to 0.
// NA years are skipped (see build_ecdf_buf). Outputs are doubles.
//
// Moved from mco-snodas src/ecdf.cpp (62fe19a). Changes: the exports are
// internal (cpp_*), every function takes `min_valid` (see ecdf_impl.h),
// cpp_ecdf_q takes probabilities as a 1 x m (shared) or n x m (row-specific)
// matrix, and ecdf_stats() is composed in R (dist_mean() + dist_quantile()).


// [[Rcpp::export]]
NumericVector cpp_ecdf_p(NumericMatrix raw, NumericVector obs, IntegerVector min_valid) {
  // raw: n_pix x n_time reference record; obs: n_pix observed values.

  int n_pix  = raw.nrow();
  int n_time = raw.ncol();

  NumericVector out(n_pix);
  std::vector<double> sorted_pos;
  sorted_pos.reserve(n_time);

  for (int i = 0; i < n_pix; i++) {
    double v = obs[i];
    if (ISNAN(v)) { out[i] = NA_REAL; continue; }

    int n_zero, n_valid;
    if (!build_ecdf_buf(raw, i, min_valid_at(min_valid, i), sorted_pos, n_zero, n_valid)) { out[i] = NA_REAL; continue; }

    out[i] = ecdf_p_pixel(sorted_pos, n_zero, n_valid, v);
  }

  return out;
}


// [[Rcpp::export]]
NumericMatrix cpp_ecdf_q(NumericMatrix raw, NumericMatrix probs, IntegerVector min_valid) {
  // raw: n_pix x n_time reference record; probs: requested probabilities,
  // 1 x n_probs (every row) or n_pix x n_probs.

  int n_pix   = raw.nrow();
  int n_time  = raw.ncol();
  int n_probs = probs.ncol();

  NumericMatrix out(n_pix, n_probs);
  std::vector<double> sorted_pos;
  sorted_pos.reserve(n_time);

  for (int i = 0; i < n_pix; i++) {
    int n_zero, n_valid;
    if (!build_ecdf_buf(raw, i, min_valid_at(min_valid, i), sorted_pos, n_zero, n_valid)) {
      for (int k = 0; k < n_probs; k++) out(i, k) = NA_REAL;
      continue;
    }
    for (int k = 0; k < n_probs; k++)
      out(i, k) = ecdf_q_pixel(sorted_pos, n_zero, n_valid, prob_at(probs, i, k));
  }

  return out;
}


// [[Rcpp::export]]
NumericVector cpp_ecdf_mean(NumericMatrix raw, IntegerVector min_valid) {
  // Mean of the valid years, zeros included; NA if too few are valid.

  int n_pix  = raw.nrow();
  int n_time = raw.ncol();
  NumericVector out(n_pix);

  for (int i = 0; i < n_pix; i++) {
    double sum  = 0.0;
    int n_valid = 0;
    for (int j = 0; j < n_time; j++) {
      double h = raw(i, j);
      if (ISNAN(h)) continue;
      sum += h;
      n_valid++;
    }
    out[i] = n_valid < min_valid_at(min_valid, i) ? NA_REAL : sum / n_valid;
  }

  return out;
}
