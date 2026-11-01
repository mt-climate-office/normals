# Today's function names (mco-snodas src/ and R/calc_snodas_zig_spline_normals.R)
# as 2.0 calls (API.md 6.1), so the ported scripts/test_distributions.R reads
# as it did. min_valid is today's rule, min(10, record length) (D3).

legacy_mv <- function(raw) min_valid_years(ncol(raw), cap = 10)
legacy_zig <- family_zig()
legacy_emp <- family_empirical()
as_mat <- function(x) if (is.matrix(x)) x else matrix(x, nrow = 1L)

ecdf_p <- function(raw, obs) {
  raw <- as_mat(raw)
  dist_cdf(normals_dist(legacy_emp, NULL, raw, legacy_mv(raw), NULL), obs)
}
# Probabilities outside [0, 1] have no quantile_columns() name, so the
# internal C++ is called directly for them (today's ecdf_q accepted them).
ecdf_q <- function(raw, probs) {
  raw <- as_mat(raw)
  if (any(probs < 0 | probs > 1))
    return(normals:::cpp_ecdf_q(normals:::as_double_matrix(raw), matrix(probs, 1L), legacy_mv(raw)))
  unname(dist_quantile(normals_dist(legacy_emp, NULL, raw, legacy_mv(raw), NULL), probs))
}
ecdf_mean <- function(raw) {
  raw <- as_mat(raw)
  dist_mean(normals_dist(legacy_emp, NULL, raw, legacy_mv(raw), NULL))
}
ecdf_stats <- function(raw) {
  raw <- as_mat(raw)
  d <- normals_dist(legacy_emp, NULL, raw, legacy_mv(raw), NULL)
  list(mean = dist_mean(d), quantiles = dist_quantile(d, NORMALS_QUANTILES))
}
zig_lmom <- function(mat) day_stats(as_mat(mat), legacy_zig)
zig_params_from_lmom <- function(lm, n_time) {
  unname(fit_day(lm, legacy_zig, min_valid_years(n_time, cap = 10)))
}
zig_fit <- function(mat) {
  mat <- as_mat(mat)
  unname(fit_day(day_stats(mat, legacy_zig), legacy_zig, legacy_mv(mat)))
}
zig_shape_from_lcv <- function(lcv) normals:::cpp_zig_shape_from_lcv(lcv)
zig_p <- function(params, obs, raw) {
  raw <- as_mat(raw)
  dist_cdf(normals_dist(legacy_zig, params, raw, legacy_mv(raw), "empirical"), obs)
}
zig_q <- function(params, probs, raw) {
  raw <- as_mat(raw)
  unname(dist_quantile(normals_dist(legacy_zig, params, raw, legacy_mv(raw), "empirical"), probs))
}
zig_mean <- function(params, raw) {
  raw <- as_mat(raw)
  dist_mean(normals_dist(legacy_zig, params, raw, legacy_mv(raw), "empirical"))
}
ZIG_SMOOTH_LCV_RANGE <- legacy_zig$control$lcv_range
zig_spline_basis <- function(k = legacy_zig$control$k, at = 1:365) smooth_basis(at, k)
zig_spline_suffstats <- function(mat) zig_lmom(mat)
# Today's zig_spline_fit_stats(stats, basis): the fit at the basis's days 1:365.
zig_spline_fit_stats <- function(stats, basis) smooth_fit(stats, seq_len(nrow(basis$X)), legacy_zig)
zig_spline_params <- function(coefs, basis, status = NULL) {
  X <- if (is.list(basis)) basis$X else basis
  normals:::zig_smooth_eval(coefs, X, status, legacy_zig$control)
}
