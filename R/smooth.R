#' Smoothing through the year
#'
#' * `smooth_basis()`: the cyclic cubic P-spline basis (mgcv `bs = "cp"`,
#'   knots 0.5 - 365.5, no constraint absorbed) at cycle positions `at`, and
#'   its penalty: `list(X, S, k)`.
#' * `stack_stats()`: a list of `day_stats()` matrices, one per day, as the
#'   stat-major list `smooth_fit()` takes (one n x n_days matrix per
#'   statistic).
#' * `smooth_fit()`: fits the family's smoother to every row. Rows are split
#'   into `cores` contiguous chunks and fitted in forked workers; the result
#'   does not depend on `cores`. Returns a `normals_smooth`: `coef`, `edf`,
#'   `iter`, `status`, `family`, `at`, `k`, `lambda`.
#' * `smooth_eval()`: the parameters at cycle positions `at`, a list (by
#'   parameter) of n x `length(at)` matrices.
#' * `smooth_params()`: the n x p parameter matrix for position `j` of a
#'   `smooth_eval()` result.
#' * `fit_status_label()`: status codes as spec section 8 `fit_status` labels.
#' * `smooth_tables()`: today's coefficient and diagnostics tables (layout
#'   v0), `list(coef, diagnostics)`.
#'
#' @param at Cycle positions ([cycle_day()]).
#' @param k Number of basis functions.
#' @param stats_by_day List of `day_stats()` matrices, one per day.
#' @param stats Stat-major daily statistics (see `stack_stats()`).
#' @param family A `normals_family` with a `smooth` slot.
#' @param cores Forked workers.
#' @param fit A `normals_smooth`.
#' @param evaluated A `smooth_eval()` result.
#' @param j Position index.
#' @param status Integer status codes.
#' @param subjects Data frame of subject key columns, one row per fit row.
#' @name smoothing
NULL

#' @rdname smoothing
#' @export
smooth_basis <- function(at, k) not_implemented()

#' @rdname smoothing
#' @export
stack_stats <- function(stats_by_day) not_implemented()

#' @rdname smoothing
#' @export
smooth_fit <- function(stats, at, family, cores = 1L) not_implemented()

#' @rdname smoothing
#' @export
smooth_eval <- function(fit, at) not_implemented()

#' @rdname smoothing
#' @export
smooth_params <- function(evaluated, j) not_implemented()

#' @rdname smoothing
#' @export
fit_status_label <- function(status, family) not_implemented()

#' @rdname smoothing
#' @export
smooth_tables <- function(fit, subjects) not_implemented()
