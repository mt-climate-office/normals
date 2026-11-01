#include "ecdf_impl.h"

// Moved from mco-snodas src/zig.cpp (62fe19a). Changes:
//  - the exports are internal (cpp_*);
//  - the L-moments' minimum of 4 positive values is the `min_pos` argument
//    (zig_control()$min_pos_fit);
//  - the minimum number of valid years is the `min_valid` argument (see
//    ecdf_impl.h), where it was min(10, n_time);
//  - cpp_zig_q takes probabilities as a 1 x m (shared) or n x m
//    (row-specific) matrix;
//  - cpp_zig_p, cpp_zig_q and cpp_zig_mean take `fallback`: TRUE is today's
//    ECDF fallback, FALSE gives NA on the rows that would fall back.
// Arithmetic and evaluation order are unchanged.

// ── L-moments ───────────────────────────────────────────────────────────────
// The ZIG chain: a reference record → L-moments of its positive values →
// gamma shape and scale. zig_fit() runs it per day; the ZIG-smooth normals
// (src/zig_spline.cpp) smooth the same per-day L-moments over the year.

struct LmomRow {
  int    n_valid;   // valid (non-NA) years
  int    n_zero;    // zero years
  double l1;        // mean of the positive values
  double lcv;       // L-CV, l2 / l1
};

// L-moments of row i's positive values (buf is scratch space). l1 and lcv are
// NA with fewer than min_pos (4) positive values or a degenerate sample (l1 or
// l2 <= 0).
static void lmom_row(const NumericMatrix& mat, int i, int min_pos, std::vector<double>& buf, LmomRow& out) {
  const int n_time = mat.ncol();
  buf.clear();
  out.n_valid = 0;
  out.n_zero  = 0;
  out.l1      = NA_REAL;
  out.lcv     = NA_REAL;
  for (int j = 0; j < n_time; j++) {
    const double v = mat(i, j);
    if (ISNAN(v)) continue;
    out.n_valid++;
    if (v == 0.0) out.n_zero++;
    else if (v > 0.0) buf.push_back(v);
  }
  const int n_pos = buf.size();
  if (n_pos < min_pos) return;

  std::sort(buf.begin(), buf.end());

  // Unbiased sample L-moments from probability-weighted moments (Hosking 1996):
  //   b0 = mean(x), b1 = mean(x * (rank - 1) / (n - 1)), l1 = b0, l2 = 2 b1 - b0
  double b0 = 0.0, b1 = 0.0;
  for (int j = 0; j < n_pos; j++) {
    b0 += buf[j];
    b1 += buf[j] * (double)j / (n_pos - 1);
  }
  b0 /= n_pos;
  b1 /= n_pos;

  const double l1 = b0;
  const double l2 = 2.0 * b1 - b0;
  if (l2 <= 0.0 || l1 <= 0.0) return;
  out.l1  = l1;
  out.lcv = l2 / l1;
}

// Gamma shape from the L-CV: Hosking's rational approximation (PELGAM,
// LMOMENTS Fortran package).
inline double shape_from_lcv(double cv) {
  const double A1 = -0.3080;
  const double A2 = -0.05812;
  const double A3 =  0.01765;
  const double B1 =  0.7213;
  const double B2 = -0.5947;
  const double B3 = -2.1817;
  const double B4 =  1.2113;
  if (cv >= 0.5) {
    const double t = 1.0 - cv;
    return t * (B1 + t * B2) / (1.0 + t * (B3 + t * B4));
  }
  const double t = M_PI * cv * cv;
  return (1.0 + A1 * t) / (t * (1.0 + t * (A2 + t * A3)));
}

// One row's (shape, scale, p_zero) from its L-moments; all NA with fewer than
// min_valid valid years, shape and scale NA without L-moments.
inline void params_from_lmom(const LmomRow& lm, int min_valid, double& shape, double& scale, double& p_zero) {
  shape = scale = p_zero = NA_REAL;
  if (lm.n_valid < min_valid) return;
  p_zero = (double)lm.n_zero / lm.n_valid;
  if (ISNAN(lm.l1)) return;
  shape = shape_from_lcv(lm.lcv);
  scale = lm.l1 / shape;   // as in R's pgamma(scale =)
}

// [[Rcpp::export]]
NumericMatrix cpp_zig_lmom(NumericMatrix mat, int min_pos) {
  // Per-row L-moment summary of a reference record: the first link of the
  // ZIG chain, and the per-day statistics of the ZIG-smooth fit.
  //
  // mat: n_pix x n_time reference record (NA = missing year).
  // Returns n_pix x 4: n (valid years), z (zero years), l1 (mean of the
  // positive values) and lcv (their L-CV, l2 / l1). l1 and lcv are NA with
  // fewer than 4 positive values or a degenerate sample; n and z are counted
  // whatever the record length.
  const int n_pix = mat.nrow();
  NumericMatrix out(n_pix, 4);
  std::vector<double> buf;
  buf.reserve(mat.ncol());
  LmomRow lm;
  for (int i = 0; i < n_pix; i++) {
    lmom_row(mat, i, min_pos, buf, lm);
    out(i, 0) = lm.n_valid;
    out(i, 1) = lm.n_zero;
    out(i, 2) = lm.l1;
    out(i, 3) = lm.lcv;
  }
  colnames(out) = CharacterVector::create("n", "z", "l1", "lcv");
  return out;
}

// [[Rcpp::export]]
NumericMatrix cpp_zig_params_from_lmom(NumericMatrix lmom, IntegerVector min_valid) {
  // ZIG parameters from zig_lmom() output, as zig_fit() computes them.
  //
  // lmom:      n_pix x 4 (n, z, l1, lcv) from zig_lmom().
  // min_valid: minimum valid years, length 1 or n_pix.
  // Returns n_pix x 3 (shape, scale, p_zero).
  const int n_pix = lmom.nrow();
  NumericMatrix out(n_pix, 3);
  LmomRow lm;
  for (int i = 0; i < n_pix; i++) {
    lm.n_valid = lmom(i, 0);
    lm.n_zero  = lmom(i, 1);
    lm.l1      = lmom(i, 2);
    lm.lcv     = lmom(i, 3);
    params_from_lmom(lm, min_valid_at(min_valid, i), out(i, 0), out(i, 1), out(i, 2));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector cpp_zig_shape_from_lcv(NumericVector lcv) {
  // Gamma shape from L-CV (Hosking's approximation), elementwise; NA stays NA.
  NumericVector out(lcv.size());
  for (R_xlen_t i = 0; i < lcv.size(); i++)
    out[i] = ISNAN(lcv[i]) ? NA_REAL : shape_from_lcv(lcv[i]);
  return out;
}

// [[Rcpp::export]]
NumericMatrix cpp_zig_fit(NumericMatrix mat, int min_pos, IntegerVector min_valid) {
  // Zero-inflated gamma fit per row: p_zero is the fraction of zero years;
  // shape and scale are fit to the positive values by L-moments
  // (zig_lmom() → zig_params_from_lmom()).
  //
  // mat: n_pix x n_time reference record.
  // Returns n_pix x 3 (shape, scale, p_zero). NA years are skipped; all three
  // are NA with too few valid years (see min_valid_years). shape and scale are
  // NA with fewer than 4 positive values or a degenerate fit (l1 or l2 <= 0).
  // (Kept internal to test that it equals the chain.)
  const int n_pix  = mat.nrow();
  const int n_time = mat.ncol();
  NumericMatrix out(n_pix, 3);
  std::vector<double> buf;
  buf.reserve(n_time);
  LmomRow lm;
  for (int i = 0; i < n_pix; i++) {
    lmom_row(mat, i, min_pos, buf, lm);
    params_from_lmom(lm, min_valid_at(min_valid, i), out(i, 0), out(i, 1), out(i, 2));
  }
  return out;
}

// Zero-inflated gamma distribution:
//   P(X = 0)  = p_zero
//   P(X <= x) = p_zero + (1 - p_zero) * pgamma(x, shape, scale),   x > 0
//
// zig_p maps a zero observation to the mid-probability p_zero / 2, so a zero
// in an all-zero record is the 50th percentile, not the 100th. zig_q returns
// 0 for p <= p_zero, so zig_q(zig_p(x)) = x for x = 0 and, on the gamma
// path, for x > 0.
//
// Where shape or scale is NA (fewer than 4 positive values or a degenerate
// fit), both fall back to the zero-inflated type-7 ECDF of raw (ecdf_impl.h),
// skipping NA years.
// NA p_zero or obs gives NA, as does p_adj >= 1 in zig_q (so p = 1 is NA).
// Outputs are doubles.

// [[Rcpp::export]]
NumericVector cpp_zig_p(NumericMatrix params,
                        NumericVector obs,
                        NumericMatrix raw,
                        IntegerVector min_valid,
                        bool fallback) {
  // params:    n_pix x 3 (shape, scale, p_zero), as from zig_fit().
  // obs:       n_pix observed values.
  // raw:       n_pix x n_time reference record, for the ECDF fallback.
  // min_valid: minimum valid years of the fallback, length 1 or n_pix.
  // fallback:  FALSE gives NA where the ECDF fallback would be used.

  int n_pix  = params.nrow();
  int n_time = raw.ncol();

  NumericVector out(n_pix);
  std::vector<double> sorted_pos;
  sorted_pos.reserve(n_time);

  for (int i = 0; i < n_pix; i++) {
    double shape  = params(i, 0);
    double scale  = params(i, 1);
    double p_zero = params(i, 2);
    double v      = obs[i];

    if (ISNAN(p_zero) || ISNAN(v)) { out[i] = NA_REAL; continue; }

    // Zero observation: mid-probability.
    if (v == 0.0) { out[i] = p_zero / 2.0; continue; }

    if (!ISNAN(shape) && !ISNAN(scale)) {
      out[i] = p_zero + (1.0 - p_zero) * R::pgamma(v, shape, scale, 1, 0);
      continue;
    }

    // ECDF fallback.
    int n_zero, n_valid;
    if (!fallback || !build_ecdf_buf(raw, i, min_valid_at(min_valid, i), sorted_pos, n_zero, n_valid)) {
      out[i] = NA_REAL;
    } else {
      out[i] = ecdf_p_pixel(sorted_pos, n_zero, n_valid, v);
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericMatrix cpp_zig_q(NumericMatrix params,
                        NumericMatrix probs,
                        NumericMatrix raw,
                        IntegerVector min_valid,
                        bool fallback) {
  // params:    n_pix x 3 (shape, scale, p_zero), as from zig_fit().
  // probs:     requested probabilities, 1 x n_probs (every row) or
  //            n_pix x n_probs.
  // raw:       n_pix x n_time reference record, for the ECDF fallback.
  // min_valid: minimum valid years of the fallback, length 1 or n_pix.
  // fallback:  FALSE gives NA where the ECDF fallback would be used.
  // Returns n_pix x n_probs quantiles (mm).

  int n_pix   = params.nrow();
  int n_probs = probs.ncol();
  int n_time  = raw.ncol();

  NumericMatrix out(n_pix, n_probs);

  std::vector<double> sorted_pos;
  sorted_pos.reserve(n_time);

  for (int i = 0; i < n_pix; i++) {
    double shape  = params(i, 0);
    double scale  = params(i, 1);
    double p_zero = params(i, 2);

    if (ISNAN(p_zero)) {
      for (int k = 0; k < n_probs; k++) out(i, k) = NA_REAL;
      continue;
    }

    bool use_empirical = ISNAN(shape) || ISNAN(scale);

    // ECDF fallback buffer, built once per row.
    bool too_few    = false;
    int  n_zero_fb  = 0;
    int  n_valid_fb = 0;
    if (use_empirical) {
      too_few = !fallback || !build_ecdf_buf(raw, i, min_valid_at(min_valid, i), sorted_pos, n_zero_fb, n_valid_fb);
    }

    for (int k = 0; k < n_probs; k++) {
      double p     = prob_at(probs, i, k);
      double p_adj = (p - p_zero) / (1.0 - p_zero);

      if (p_adj <= 0.0) { out(i, k) = 0.0;     continue; }
      if (p_adj >= 1.0) { out(i, k) = NA_REAL;  continue; }

      if (!use_empirical) {
        out(i, k) = R::qgamma(p_adj, shape, scale, 1, 0);
        continue;
      }

      // ECDF fallback; p <= p_zero was already mapped to 0 above.
      out(i, k) = too_few ? NA_REAL : ecdf_q_pixel(sorted_pos, n_zero_fb, n_valid_fb, p);
    }
  }

  return out;
}

// Mean of the zero-inflated gamma: E[X] = (1 - p_zero) * shape * scale.
// NA p_zero gives NA; NA shape or scale falls back to the mean of the valid
// years of raw, zeros included.

// [[Rcpp::export]]
NumericVector cpp_zig_mean(NumericMatrix params, NumericMatrix raw,
                           IntegerVector min_valid, bool fallback) {
  // params:    n_pix x 3 (shape, scale, p_zero), as from zig_fit().
  // raw:       n_pix x n_time reference record, for the empirical fallback.
  // min_valid: minimum valid years of the fallback, length 1 or n_pix.
  // fallback:  FALSE gives NA where the empirical fallback would be used.

  int n_pix  = params.nrow();
  int n_time = raw.ncol();
  NumericVector out(n_pix);

  for (int i = 0; i < n_pix; i++) {
    double shape  = params(i, 0);
    double scale  = params(i, 1);
    double p_zero = params(i, 2);

    if (ISNAN(p_zero)) {
      out[i] = NA_REAL;
      continue;
    }

    if (!ISNAN(shape) && !ISNAN(scale)) {
      out[i] = (1.0 - p_zero) * shape * scale;
      continue;
    }

    if (!fallback) {
      out[i] = NA_REAL;
      continue;
    }

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
