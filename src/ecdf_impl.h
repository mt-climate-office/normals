#pragma once
#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Zero-inflated type-7 ECDF helpers, shared by ecdf.cpp and zig.cpp (ECDF
// fallback). No Rcpp exports here; they belong in the .cpp files.
//
// Moved from mco-snodas src/ecdf_impl.h (62fe19a). The one change: the
// minimum number of valid years is the caller's `min_valid` argument (an
// integer, length 1 or one per row), where it was min(10, n_time) here.


// NA years are skipped: SNODAS's valid-data mask grew over the record (mostly
// water bodies), so some pixels and elevation bands are NA in early years
// only. A row needs min_valid valid years, else NA.

// Row i's minimum: min_valid is length 1 (every row) or one per row.
inline int min_valid_at(const IntegerVector& min_valid, int i) {
  return min_valid.size() == 1 ? min_valid[0] : min_valid[i];
}

// Row i's probability k: probs is 1 x m (every row) or n x m (row-specific).
inline double prob_at(const NumericMatrix& probs, int i, int k) {
  return probs.nrow() == 1 ? probs(0, k) : probs(i, k);
}

// For row i: sorted positive values, zero count and valid-year count.
// Returns false (caller outputs NA) if too few years are valid. n_valid, not
// n_time, is the sample size of the distribution.
static bool build_ecdf_buf(NumericMatrix raw, int i, int min_valid,
                            std::vector<double>& sorted_pos,
                            int& n_zero, int& n_valid) {
  int n_time = raw.ncol();
  sorted_pos.clear();
  n_zero  = 0;
  n_valid = 0;
  for (int j = 0; j < n_time; j++) {
    double h = raw(i, j);
    if (ISNAN(h)) continue;
    n_valid++;
    if (h == 0.0) n_zero++;
    else          sorted_pos.push_back(h);
  }
  if (n_valid < min_valid) return false;
  std::sort(sorted_pos.begin(), sorted_pos.end());
  return true;
}


// Zero-inflated type-7 CDF for one row; n_time is the valid-year count.
// v == 0 returns the mid-probability p_zero / 2; v below the smallest positive
// returns p_zero, above the largest returns 1; an all-zero record returns 1
// for any positive v.
inline double ecdf_p_pixel(const std::vector<double>& sorted_pos,
                            int n_zero, int n_time, double v) {
  int    n_pos  = (int)sorted_pos.size();
  double p_zero = (double)n_zero / n_time;

  if (v == 0.0)  return p_zero / 2.0;
  if (n_pos == 0) return 1.0;   // all-zero record: positive obs is unprecedented

  int lo = (int)(std::upper_bound(sorted_pos.begin(), sorted_pos.end(), v)
                 - sorted_pos.begin()) - 1;
  double p_adj;
  if (lo < 0) {
    p_adj = 0.0;
  } else if (lo >= n_pos - 1) {
    p_adj = 1.0;
  } else {
    double span = sorted_pos[lo + 1] - sorted_pos[lo];
    double frac = (span > 0.0) ? (v - sorted_pos[lo]) / span : 0.0;
    p_adj = (lo + frac) / (n_pos - 1);
  }
  return p_zero + (1.0 - p_zero) * p_adj;
}


// Zero-inflated type-7 quantile for one row; n_time is the valid-year count.
// p < p_zero returns 0. The strict < matches R's quantile(type = 7): p = 0 on
// an all-positive row returns the smallest value, and p = p_zero returns the
// smallest positive value. zig_q uses <= and returns 0 at p = p_zero.
// p > 1 returns NA; an all-zero record returns 0 for p <= 1.
inline double ecdf_q_pixel(const std::vector<double>& sorted_pos,
                            int n_zero, int n_time, double p) {
  int    n_pos  = (int)sorted_pos.size();
  double p_zero = (double)n_zero / n_time;

  if (n_pos == 0) return (p <= 1.0) ? 0.0 : NA_REAL;

  double p_adj = (p - p_zero) / (1.0 - p_zero);
  if (p_adj < 0.0)  return 0.0;
  if (p_adj > 1.0)  return NA_REAL;

  double rank_float = p_adj * (n_pos - 1);
  int lo = (int)rank_float;
  lo = std::max(0, std::min(lo, n_pos - 1));
  int hi = std::min(lo + 1, n_pos - 1);
  double frac = rank_float - (double)lo;
  frac = std::max(0.0, std::min(frac, 1.0));

  return sorted_pos[lo] + frac * (sorted_pos[hi] - sorted_pos[lo]);
}
